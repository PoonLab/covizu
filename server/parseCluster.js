// regular expression to remove redundant sequence name components
const pat = /^hCoV-19\/(.+\/.+)\/20[0-9]{2}$/gi;
const { unique, mode, tabulate, merge_tables, utcDate } = require('./utils')
const {$DATA_FOLDER} = require("../config/config")
const dbstats = require(`../${$DATA_FOLDER}/dbstats.json`)
const d3 = require('../js/d3')
require("../globalVariables")

/**
 * Parse nodes that belong to the same variant.
 * A variant is a collection of genomes that are indistinguishable with respect to
 * (1) differences from the reference or (2) placement in the phylogenetic tree.
 * Visually, the variant is represented by a horizontal line segment spanning the
 * sample collection dates.
 * The samples that comprise a variant are gathered by collection date into "points"
 * along the horizontal line.  If there are no samples, then the variant is
 * "unsampled" and spans the entire width of the beadplot.
 *
 * @param {Object} variant:  associative list member of cluster.nodes
 * @param {number} y:  vertical position
 * @param {number} cidx:  cluster index for labeling points
 * @param {string} accn:  accession of baseline sample of variant
 * @param {Date} mindate:  used only for unsampled variants
 * @param {Date} maxdate:  used only for unsampled variants
 * @returns {{variants: [], points: []}}
 */
const parse_variant = (variant, y, cidx, accn, mindate, maxdate) => {
  var vdata, pdata = [];

  if (variant.length === 0) {
    // handle unsampled internal node
    if (mindate === null || maxdate === null) {
      // this would happen if a cluster comprised only one unsampled node - should be impossible!
      console.log("Error in parse_variants(): cannot draw unsampled variant without min and max dates");
    }
    vdata = {
      'accession': accn,
      'label': accn,
      'x1': utcDate(mindate),  // cluster min date
      'x2': utcDate(maxdate),  // cluster max date
      'y1': y,
      'y2': y,
      'count': 0,
      'country': null,
      'region': null,
      'numBeads': 0,
      'parent': null,
      'dist': 0,
      'unsampled': true
    };
  }
  else {
    // parse samples within variant, i.e., "beads"
    // [0: coldate, 1: accession, 2: location, 3: name]
    var label = variant[0][3].replace(pat, "$1"),
        coldates = variant.map(x => x[0]),
        isodate, samples, regions;

    coldates.sort();
    //Retrieving countries from variants?
    var location = variant.map(x => x[2].split(' / ')),
        region = location.map(x => x[0]),
        country = location.map(x => x[1]),
        isodates = unique(coldates);

    // update country to region map
    for (let i=0; i < country.length; i++) {
      let this_country = country[i],
          this_region = region[i];
      if (global.region_map[this_country] === undefined) {
        global.region_map[this_country] = this_region;
      }
    }

    vdata = {
      'accession': accn,
      'label': label,
      'x1': utcDate(coldates[0]),  // min date
      'x2': utcDate(coldates[coldates.length-1]),  // max date
      'y1': y,
      'y2': y,
      'count': coldates.length,
      'country': tabulate(country),
      'region': tabulate(region),
      'numBeads': isodates.length,
      'parent': null,
      'dist': 0,
      'unsampled': false
    };

    for (var i=0; i<isodates.length; i++) {
      isodate = isodates[i];
      samples = variant.filter(x => x[0] === isodate);
      location = samples.map(x => x[2].split(' / '));
      region = location.map(x => x[0]);
      country = location.map(x => x[1]);

      pdata.push({
        cidx,
        'variant': accn,
        'x': utcDate(isodate),
        'y': y,
        'count': samples.length,
        'accessions': samples.map(x => x[1]),
        'labels': samples.map(x => x[2].replace(pat, "$1")),
        'region1': mode(region),
        'region': tabulate(region),
        'country': tabulate(country),
        'parent': null,
        'dist': 0
      })
    }
  }

  return {'variant': vdata, 'points': pdata};
}


/**
 *
 * @param cluster
 * @param variants
 * @param points
 * @returns {[]}
 */
const parse_edgelist = (cluster, variants, points) => {
  // map earliest collection date of child node to vertical edges
  let edge, parent, child, dist, support,
      edgelist = [];

  // generate maps of variants and points keyed by accession
  let lookup_variant = {};
  variants.forEach(function(row) {
    lookup_variant[row.accession] = row;
  });

  let lookup_points = {},
      index_points = {};  // index by y-coordinate

  points.forEach(function(pt) {
    if (lookup_points[pt.variant] === undefined) {
      lookup_points[pt.variant] = [];
    }
    lookup_points[pt.variant].push(pt);
    if (index_points[pt.y] === undefined) {
      index_points[pt.y] = [];
    }
    index_points[pt.y].push(pt);
  })

  for (var e = 0; e < cluster.edges.length; e++) {
    edge = cluster.edges[e];
    //parent = variants.filter(x => x.accession === edge[0])[0];
    parent = lookup_variant[edge[0]];
    //child = variants.filter(x => x.accession === edge[1])[0];
    child = lookup_variant[edge[1]];

    if (parent === undefined || child === undefined) {
      // TODO: handle edge to unsampled node
      continue;
    }

    dist = parseFloat(edge[2]);
    if (edge[3] === null) {
      support = undefined;
    } else {
      support = parseFloat(edge[3]);
    }

    edgelist.push({
      'y1': parent.y1,
      'y2': child.y1,
      'x1': child.x1,  // vertical line segment
      'x2': child.x1,
      'parent': parent.label,
      'child': child.label,
      'dist': dist,
      'support': support
    });

    child.parent = parent.label;
    child.dist = dist;

    // Assign the parent and genomic distance of each point
    if (index_points[child.y1] !== undefined) {
      for (let pt of index_points[child.y1]) {
        pt.parent = parent.label;
        pt.dist = dist;
      }
    }
    /*
    for (let c = 0; c < points.length; c++) {
      if (points[c].y === child.y1) {
        points[c].parent = parent.label;
        points[c].dist = dist;
      }
    }
    */

    // update variant time range
    if (parent.x1 > child.x1) {
      parent.x1 = child.x1;
    }
    if (parent.x2 < child.x1) {
      parent.x2 = child.x1;
    }
  }
  return edgelist;
}


/**
 * Parse node and edge data from clusters JSON to a format that is
 * easier to map to SVG.
 * @param {Object} clusters:
 */
const parse_clusters = (clusters) => {
  var cluster, variant, coldates, regions, labels,
      accn, mindate, maxdate,
      result, vdata, pdata,
      variants,  // horizontal line segment data + labels
      edgelist,  // vertical line segment data
      points,  // the "beads"
      beaddata = [];  // return value

  for (const cidx in clusters) {
    cluster = clusters[cidx];

    variants = [];
    points = [];
    edgelist = [];

    // deal with edge case of cluster with only one variant, no edges
    if (Object.keys(cluster["nodes"]).length === 1) {
      variant = Object.values(cluster.nodes)[0];
      accn = Object.keys(cluster.nodes)[0];
      result = parse_variant(variant, 0, cidx, accn, null, null);
      vdata = result['variant'];
      variants.push(vdata);

      pdata = result['points'];
      points = points.concat(pdata);
    }
    else {
      // de-convolute edge list to get node list in preorder
      var nodelist = unique(cluster.edges.map(x => x.slice(0, 2)).flat());

      // date range of cluster
      coldates = nodelist.map(a => cluster.nodes[a].map(x => x[0])).flat();
      coldates.sort()
      mindate = coldates[0];
      maxdate = coldates[coldates.length - 1];

      // extract the date range for each variant in cluster
      var y = 1;
      for (const accn of nodelist) {
        // extract collection dates for all samples of this variant
        variant = cluster.nodes[accn];
        result = parse_variant(variant, y, cidx, accn, mindate, maxdate);
        variants.push(result['variant']);
        points = points.concat(result['points']);
        y++;
      }

      edgelist = parse_edgelist(cluster, variants, points);
    }

    beaddata.push({
      'variants': variants,
      'edgelist': edgelist,
      'points': points
    });

    // calculate consensus region for cluster
    // collect all region Arrays for all samples, all variants
    cluster['allregions'] = merge_tables(points.map(x => x.region));
    let max_freq = 0, max_region = '';
    for (let row of Object.entries(cluster['allregions'])) {
      if (row[1] > max_freq) {
        max_freq = row[1];
        max_region = row[0];
      }
    }
    cluster['region'] = max_region;  // most common region

    // concatenate all sample labels within cluster for searching
    labels = points.map(x => x.labels).flat();

    // decompose labels and only keep unique substrings
    let uniq = new Set(labels.map(x => x.split('/')).flat());
    cluster['searchtext'] = Array.from(uniq).join();
    cluster['label1'] = labels[0];

    // collect all countries
    cluster['country'] = merge_tables(variants.map(x => x.country));
  }
  return beaddata;
}


/**
 * Converts Date to x-coordinate of the tree scale
 * 
 * @param {Date} coldate
 * @returns {float} x-coordinate value 
 */
 function date_to_xaxis(reftip, coldate) {
  var numDays = d3.timeDay.count(reftip.first_date, coldate)
  return (numDays/365.25) + reftip.x;
}


/**
 * Map cluster information to tips of the tree.
 * @param {Array} df: data frame extracted from time-scaled tree
 * @param {Array} clusters: data from clusters JSON
 * @returns {Array} subset of data frame annotated with cluster data
 */
const map_clusters_to_tips = (df, clusters) => {
  let recombinants = [];
  // extract accession numbers from phylogeny data frame
  var tips = df.filter(x => x.children.length===0),
      tip_labels = tips.map(x => x.thisLabel),  // accessions
      tip_stats;

  for (const cidx in clusters) {
    var cluster = clusters[cidx];
    if (cluster["nodes"].length === 1) {
      continue
    }

    // find variant in cluster that matches a tip label
    var labels = Object.keys(cluster["nodes"]),
        root = tip_labels.filter(value => value === cluster['lineage'])[0];
   
    if (dbstats["lineages"][cluster["lineage"]]["raw_lineage"].startsWith("X")) {
      recombinants.push(cidx)
      continue;
    }

    if (root === undefined) {
        console.log("Failed to match cluster of index ", cidx, " to a tip in the tree");
        continue;
    }

    var root_idx = tip_labels.indexOf(root),  // row index in data frame
        root_xcoord = tips[root_idx].x;  // left side of cluster starts at end of tip

    // find most recent sample collection date
    var coldates = Array(),
        label, variant;

    for (var i=0; i<labels.length; i++) {
      label = labels[i];
      variant = cluster['nodes'][label];
      coldates = coldates.concat(variant.map(x => x[0]));
    }
    coldates.sort();  // in place, ascending order

    var first_date = utcDate(coldates[0]),
        last_date = utcDate(coldates[coldates.length-1]);
    // console.log(first_date, last_date)
    
    // Calculate the mean collection date
    let date_diffs = coldates.map(x => (utcDate(x)-first_date) / (1000 * 3600 * 24)
            // d3.timeDay.count(first_date, utcDate(x))
        ),
        mean_date = Math.round(date_diffs.reduce((a, b) => a + b, 0) / date_diffs.length);

    // let date_diffs = coldates.map(x => d3.timeDay.count(first_date, utcDate(x))),
    //     mean_date = Math.round(date_diffs.reduce((a, b) => a + b, 0) / date_diffs.length);

    // augment data frame with cluster data
    tips[root_idx].cluster_idx = cidx;
    tips[root_idx].allregions = cluster.allregions;
    tips[root_idx].region = cluster.region;
    tips[root_idx].country = cluster.country;
    tips[root_idx].searchtext = cluster.searchtext;
    tips[root_idx].label1 = cluster["lineage"];
    tips[root_idx].count = coldates.length;
    tips[root_idx].varcount = cluster["sampled_variants"]; // Number for sampled variants
    tips[root_idx].sampled_varcount = labels.filter(x => x.substring(0,9) !== "unsampled").length;
    tips[root_idx].first_date = first_date;
    tips[root_idx].last_date = last_date;
    tips[root_idx].pdist = cluster.pdist;
    tips[root_idx].rdist = cluster.rdist;

    tips[root_idx].coldate = last_date;
    tips[root_idx].x1 = root_xcoord - ((last_date - first_date) / 3.154e10);
    tips[root_idx].x2 = root_xcoord;

    // map dbstats for lineage to tip
    tip_stats = dbstats["lineages"][cluster["lineage"]];
    tips[root_idx].max_ndiffs = tip_stats.max_ndiffs;
    tips[root_idx].mean_ndiffs = tip_stats.mean_ndiffs;
    tips[root_idx].nsamples = tip_stats.nsamples;
    tips[root_idx].mutations = tip_stats.mutations;
    tips[root_idx].infections = tip_stats.infections;

    // calculate residual from mean differences and mean collection date - fixes #241
    let times = coldates.map(x => utcDate(x).getTime()),
        origin = 18231,  // days between 2019-12-01 and UNIX epoch (1970-01-01)
        mean_time = times.reduce((x, y)=>x+y) / times.length / 8.64e7 - origin,
        rate = 0.0655342,  // subs per genome per day
        exp_diffs = rate * mean_time;  // expected number of differences
    tips[root_idx].residual = tip_stats.mean_ndiffs - exp_diffs;  // tip_stats.residual;

    Date.prototype.addDays = function(days) {
      var date = new Date(this.valueOf());
      date.setDate(date.getDate() + days);
      return date;
    }

    tips[root_idx].mcoldate = first_date.addDays(mean_date);
  }

  var recombinant_tips = []
  for (const cidx in recombinants) {
    var cluster = clusters[recombinants[cidx]];
    if (cluster["nodes"].length === 1) {
      continue
    }

    // find variant in cluster that matches a tip label
    var labels = Object.keys(cluster["nodes"]);
    var root_idx = cidx;  // row index in data frame

    // find most recent sample collection date
    var coldates = Array(),
        label, variant;

    for (var i=0; i<labels.length; i++) {
      label = labels[i];
      variant = cluster['nodes'][label];
      coldates = coldates.concat(variant.map(x => x[0]));
    }
    coldates.sort();  // in place, ascending order

    var first_date = utcDate(coldates[0]),
        last_date = utcDate(coldates[coldates.length-1]);
    
    var root_xcoord = date_to_xaxis(tips[0], last_date);  // left side of cluster starts at end of tip

    // Calculate the mean collection date
    let date_diffs = coldates.map(x => (utcDate(x)-first_date) / (1000 * 3600 * 24)
            // d3.timeDay.count(first_date, utcDate(x))
        ),
        mean_date = Math.round(date_diffs.reduce((a, b) => a + b, 0) / date_diffs.length);


    // calculate residual from mean differences and mean collection date - fixes #241
    let times = coldates.map(x => utcDate(x).getTime()),
        origin = 18231,  // days between 2019-12-01 and UNIX epoch (1970-01-01)
        mean_time = times.reduce((x, y)=>x+y) / times.length / 8.64e7 - origin,
        rate = 0.0655342,  // subs per genome per day
        exp_diffs = rate * mean_time;  // expected number of differences

    tip_stats = dbstats["lineages"][cluster["lineage"]];

    // recombinant_tips[root_idx].mcoldate = first_date.addDays(mean_date);
    recombinant_tips.push({
      cluster_idx: recombinants[cidx],
      allregions: cluster.allregions,
      region: cluster.region,
      country: cluster.country,
      searchtext: cluster.searchtext,
      label1: cluster["lineage"],
      count: coldates.length,
      varcount: cluster["sampled_variants"], // Number for sampled variants
      sampled_varcount: labels.filter(x => x.substring(0,9) !== "unsampled").length,
      first_date: first_date,
      last_date: last_date,
      pdist: cluster.pdist,
      rdist: cluster.rdist,
  
      coldate: last_date,
      x1: root_xcoord - ((last_date - first_date) / 3.154e10),
      x2: root_xcoord,
      y: cidx,
      
      // map dbstats for lineage to tip
      max_ndiffs: tip_stats.max_ndiffs,
      mean_ndiffs: tip_stats.mean_ndiffs,
      nsamples: tip_stats.nsamples,
      mutations: tip_stats.mutations,
      residual: tip_stats.mean_ndiffs - exp_diffs,  // tip_stats.residual;
      mcoldate: first_date.addDays(mean_date),
      infections: tip_stats.infections
    })
  }

  recombinant_tips = recombinant_tips.sort(function(a, b) {
    return a.first_date - b.first_date;
  })

  for (const cidx in recombinant_tips) {
    recombinant_tips[cidx].y = cidx
  }

  return { tips, recombinant_tips };
}


/**
 * Populate Object with accession-cluster ID as key-value pairs.
 * Note, this also provides a list (via Object.keys()) of all
 * accession numbers for autocompleting search queries.
 *
 * @param {Object} clusters:  contents of clusters JSON
 * @returns {{}}
 */
 function index_accessions(clusters) {
	var index = {};
	for (const cid in clusters) {
		var accns = Object.entries(clusters[cid].nodes)
			.map(x => x[1])
			.flat()
			.map(x => x[1]);
		for (const accn of accns) {
			index[accn] = cid;
		}
	}
	return(index); 
}


function index_lineage(clusters) {
  var index = {};
  for (const cid in clusters) {
    var accns = clusters[cid].lineage
    index[accns] = cid;
  }
  return index;
}


module.exports = {
  parse_clusters,
  map_clusters_to_tips,
  index_accessions,
  index_lineage,
};