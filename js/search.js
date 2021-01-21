// to store references to SVG objects (nodes)
var selected = [];

/**
* Prepare the search stats with initial stats
*/
function prepare_search_stats(initial_stats) {
  let stats = {
    ...initial_stats,  // ... = spread syntax, cast iterable as multiple arguments
  };

  return {
    get: () => stats,

    update: (new_stats) => {
      stats = {
        ...stats,
        ...new_stats,
      }

      return stats;
    },
  }
}

const search_stats = prepare_search_stats({
  query: undefined,
  start: undefined,
  end: undefined,
  current_point: 0,
  total_points: NaN,
  points: [],
  bead_indexer: 0,
  start_idx: []
});
 
const search_results = prepare_search_stats({
  beads: [],
  clusters: []
});

function update_search_stats(stats) {
  $('#search_stats').text(`${stats.current_point+1} of ${stats.total_points} points`);
}

/**
 * This is the main search search function, 
 * it identifies which cluster and beads contain a hit. 
 * It returns an array that contains the total number of hits 
 * and the number of hits in each cluster.
 * When the user changes the search query this should be rerun.
 *
 * @param all_bead_data: all the bead data that should be compared against the searched quary.
 * @param text_query: a text string that should be compared against the searchtext property.
 * @param start_data: a date type indicating the search start date.
 * @param end_data: a date type indicating the search end date.
 */
function main_search(all_bead_data, text_query, start_date, end_date) {
  // Flatten the json data to an array with bead data only
  flat_data = find_beads_points(all_bead_data);

  // Remove once dates are implemented
  start_date = new Date("2020-01-01");
  end_date = new Date("2021-09-30");

  //Find all the beads that are a hit
  search_hits = flat_data.filter(function(bead) {
	  temp = (bead.accessions.some(accession => accession.includes(text_query)) || 
		  bead.labels.some(label => label.includes(text_query))) && 
		  (bead.x >= start_date && bead.x <= end_date);
	  return temp;
  });

  // Move this to some where else
  var map_cidx_to_id = [], key;
  var rect = d3.selectAll('#svg-timetree > svg > rect')
	.nodes();
  for (var i = 0; i < rect.length; i++) {
	  key = d3.select(rect[i]).attr("cidx");
	  map_cidx_to_id[key] = parseInt(d3.select(rect[i]).attr("id").substring(3));
          }

  // Order the search results by cluster id, y cord, x cord 
  search_hits.sort(function(x, y) {
	  return map_cidx_to_id['cidx-' + y.cidx] - map_cidx_to_id['cidx-' + x.cidx] || 
		x.y - x.y || x.x - y.y;
  }); 
 
  // Unique identifiers for the beads that are a hit
  bead_hits = search_hits.reduce(function(map, bead, i) {
	  map[bead.accessions[0]] = i + 1; // value for the search indexing
	  return map;
  }, {});

  // Unique identifiers for the clusters that are a hit
  cluster_hits = search_hits.reduce(function(map, bead, i) {
	  map[bead.cidx] = i;
	  return map;
  }, {});

 // Update the search resutls array with the hits
 search_results.update({
	 beads: bead_hits,
	 clusters: cluster_hits
 })
}

/**
 * Highlight clusters for which the cluster function returns true
 * and highlight points in currently displayed bead plot.
 * When user switches bead plots, update highlighted points.
 *
 * @param cluster_func:  function that takes a cluster datum and returns true
 *                       when the cluster must be selected. False otherwise.
 */
function select_clusters(rects) {
  var cluster, itter;

  d3.selectAll("circle").dispatch("mouseout");  // deactivate cluster highlights

  if (d3.selectAll(".clicked").empty()) {
    // new search, jump to first matching cluster
    cluster = d3.select(rects.nodes().pop());  // last element in array
    cluster.attr("class", "clicked");
    draw_cluster_box(cluster);

    // switch to beadplot for this cluster
    beadplot(cluster.datum().cluster_idx);

    // remove all other elements from array
    itter = rects.nodes().splice(0, rects.nodes().length-1);

    d3.select("#svg-timetree")
        .selectAll("rect:not(.clicked):not(.clickedH)")
        .attr("class","not_SelectedCluster");
  }
  else {
    // FIXME: needs explanatory comment
    d3.select("#svg-timetree")
        .selectAll("rect:not(.clickedH)")
        .attr("class","not_SelectedCluster");
    itter = rects.nodes();
  }

  // FIXME: other matching clusters not getting highlighted, due to click below?
  for (const node of itter) {
    d3.select(node).attr("class", "SelectedCluster");
  }
}


/**
 * Highlight beads in current beadplot.
 * @param points_ui:  d3.Selector object
 */
function select_beads(points_ui) {
  //
  selected = points_ui.nodes();

  //
  d3.selectAll("circle:not(.selectionH)")
      .attr("class", "not_SelectedBead");
  //
  d3.select("div#svg-cluster")
      .selectAll("line")
      .attr("stroke-opacity", 0.3);
  for (const node of selected) {
    selected_obj = d3.select(node);
    create_selection(selected_obj);
  }

  // re-focus beadplot to first matching bead
  if (selected.length > 0) {
    selected[0].scrollIntoView({block: "center"});
    selected_obj = d3.select(selected[0]);
    draw_halo(selected_obj.datum());
  }
}


/**
 * Highlight clusters containing samples that match substring,
 * and highlight samples in currently displayed bead plot.
 * When user switches bead plots, update highlighted samples.
 *
 * @param {string} substr:  sub-string to search
 */
function select_beads_by_substring(substr) {
  var rects, points_ui, start_date, end_date,
    start = $("#start-date").val(),
    end = $("#end-date").val();

  if (substr === "") {
    // user submitted empty string
    clear_selection();
    return;
  }
  if (start == "") {
    start_date = new Date("2019-01-01");
  } else {
    start_date = new Date(start);  // from ISO string
  }
  if (end == "") {
    end_date = new Date();  // today
  } else {
    end_date = new Date(end);
  }

  rects = d3.selectAll("#svg-timetree > svg > rect:not(.clickedH)")
    .filter(function(d) {
      return (d.searchtext.match(substr) !== null
          //&& d.first_date <= end_date && d.last_date >= start_date
      )
    });
  main_search(beaddata, substr);
  select_clusters(rects);

  // preceding function switches view to beadplot with matching samples, if any
  points_ui = d3.selectAll("#svg-cluster > svg > g > circle")
    .filter(function(d) {
      return (d.labels.some(x => x.includes(substr))
          //&& d.x >= start_date && d.x <= end_date
      )
    });
  select_beads(points_ui);
}


/**
 * Highlight clusters containing the points.
 *
 * @param points an array of points
 */
function select_by_points(points){
  if (points.length <= 0) {
    // user submitted empty string
    clear_selection();
    return;
  }

  const labels = new Set(points.map(point => point.labels).flat());
  const cidx = new Set(points.map(point => point.cidx));
  const clusters_fn = (cluster) => cluster.cluster_idx && cidx.has(cluster.cluster_idx);
  const points_fn = (point) => point && point.labels.some(label => labels.has(label));
  select_beads_by_comparators(clusters_fn, points_fn);
}


/**
 * Highlight and jump to node in beadplot by sample accession number.
 * @param {string} accn:  accession number to search for
 * @param {Boolean} reset_clusters_tree: when false, the cluster tree remain unchanged
 */
function select_bead_by_accession(accn, reset_clusters_tree = true) {
  // switch to cluster beadplot
  var cid = accn_to_cid[accn];

  if (cid !== undefined) {

    if (reset_clusters_tree === true) {
      d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)").attr("class", "not_SelectedCluster");
    }

    var rect = d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)")
      .filter(function(d) { return(d.cluster_idx === cid); })
      .attr("class", "clicked");

    d3.select(rect.nodes().pop()).each(function(d) {
      draw_cluster_box(d, this);
    });

    beadplot(cid);

    var bead = d3.selectAll("circle").filter(function(d) {
      return d.accessions.includes(accn);
    });

    create_selection(bead);
    bead.node().scrollIntoView({block: "center"});

    // Update information panel
    const datum = bead.datum();
    gentable(datum);
    draw_region_distribution(tabulate(datum.region));
    gen_details_table(datum);

  }
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
			.map(x => x.accession);
		for (const accn of accns) {
			index[accn] = cid;
		}
	}
	return(index);
}

function as_label(search_data) {
	const [, accn] = search_data;
	return accn;
}

/**
 * Provides a source function suitable for jQuery UI's autocomplete
 *
 * @param {object} accn_to_cid: Accession numbers - cluster ids mapping
 * @returns {function}
 */
function get_autocomplete_source_fn(accn_to_cid) {
	// This is a hack to match anything that could be an acc number prefix
	const prefix = /^(E|I|EP|IS|EPI_I|EPI_IS|EPI_ISL_?|EPI_?|ISL_?)$/i;
	const MIN_RESULTS = 10;
	const normalize = (str) => str.replace(/[^a-z0-9]/gi, '').toLowerCase();
	const data = Object.keys(accn_to_cid).map(accn => [
		normalize(accn), accn
	]);

	return function({ term }, response) {
		if (!/\d/.test(term)) {
			if (prefix.test(term)) {
				response(data.slice(0, MIN_RESULTS).map(as_label));
			} else {
				response([]);
			}
		} else {
			const result = data.filter(array => array[0].indexOf(normalize(term)) > -1);
			response(result.slice(0, MIN_RESULTS).map(as_label));
		}
	}
}


// FIXME: no argument needed here
function search() {
	var query = $('#search-input').val();

  // FIX ME: Accn search returning 0 beads
  //const points = find_beads_points(beaddata)
    //.filter(point => point.labels.some(label => label.includes(query)));
  // TODO: Make select_bead_by_* use find_beads_points result
  isAccn(query) ? select_bead_by_accession(query) : select_beads_by_substring(query);
  //const stats = search_stats.update({
    //query,
    //current_point: 0,
    //total_points: points.length,
    //points: points,
    //bead_indexer: 0,
  //});
  //update_search_stats(stats);
}


function find_beads_points(beadsdata){
  return beadsdata.map(bead => bead.points).flat()
}


function search_by_dates(start, end) {
  console.log(start);
  console.log(end);
  var rects = d3.selectAll("#svg-timetree > svg > rect:not(.clickedH)")
    .filter(function(d) { return d.first_date <= end && d.last_date >= start });
  select_clusters(rects);

  // preceding function switches view to beadplot with matching samples, if any
  var points_ui = d3.selectAll("#svg-cluster > svg > g > circle")
    .filter(function(d) { return d.x >= start && d.x <= end });
  select_beads(points_ui);

  /*
  const points = find_beads_points(beaddata)
    .filter(point => point.x >= start && point.x <= end);

  select_by_points(points);

  const stats = search_stats.update({
    start,
    end,
    current_point: 0,
    total_points: points.length,
    points: points,
    bead_indexer: 0,
  });
  update_search_stats(stats);
   */
}
