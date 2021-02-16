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
  clusters_first_bead: [],
  clusters_last_bead: [],
  current_point: 0,
  total_points: 0,
  hit_ids: []
});

function update_search_stats(stats) {
  $('#search_stats').text(`${stats.current_point+1} of ${stats.total_points} points`);
}

/**
 * This is the main search search function, 
 * it identifies which cluster and beads contain a hit. 
 * It returns an array that contains the total number of hits 
 * and the number of hits in each cluster.
 *
 * @param all_bead_data: all the bead data that should be compared against the searched quary.
 * @param text_query: a text string that should be compared against the searchtext property.
 * @param start_data: a date type indicating the search start date.
 * @param end_data: a date type indicating the search end date.
 */
function main_search(all_bead_data, text_query, start_date, end_date) {
  // Flatten the json data to an array with bead data only
  flat_data = find_beads_points(all_bead_data);

  //Find all the beads that are a hit. Convert text_query to lower case and checks to see if there is a match
  search_hits = flat_data.filter(function(bead) {
	  temp = (bead.accessions.some(accession => (accession.toLowerCase()).includes(text_query.toLowerCase())) || 
		  bead.labels.some(label => (label.toLowerCase()).includes(text_query.toLowerCase()))) && 
		  (bead.x >= start_date && bead.x <= end_date);
	  return temp;
  });

  // If there are no hits, then stops the main_search
  if (search_hits.length === 0) {
    $('#error_message').text(`No matches. Please try again.`);
    return;
  }

  // Order the search results by cluster id, y cord, x cord 
  search_hits.sort(function(x, y) {
	  return map_cidx_to_id['cidx-' + y.cidx] - map_cidx_to_id['cidx-' + x.cidx] || 
		x.y - x.y || x.x - y.y;
  }); 
 
  // Unique identifiers for the beads that are a hit
  bead_hits = search_hits.reduce(function(map, bead, i) {
	  map[bead.accessions[0]] = i; // value for the search indexing
	  return map;
  }, {});

  // Dictionary to find the index of clusters that are a hit (i.e. first bead in cluster)
  var cluster_hits = [], cluster_hits_last_id = [];
  search_hits.forEach(function(bead) {
    if (cluster_hits['cidx-' + bead.cidx] === undefined)
      cluster_hits['cidx-' + bead.cidx] = bead.accessions[0];

    cluster_hits_last_id['cidx-' + bead.cidx] = bead.accessions[0];
  });

  // Stores a sorted list of the search hit cluster ids
  var hit_keys = Object.keys(cluster_hits);
  bead_id_to_accession = Object.keys(bead_hits);

  var hit_id = []
  for (var i = 0; i < hit_keys.length; i++) {
    if (map_cidx_to_id[hit_keys[i]] === undefined)
      console.log("undefined: " + hit_keys[i]);
    hit_id[i] = map_cidx_to_id[hit_keys[i]];
  }

  hit_id.sort(function(a, b) {
    return a - b;
  });

  // Keeps track of the clicked cluster
  var curr_cluster = d3.selectAll(".clicked").nodes()[0].attributes.cidx.nodeValue;
  var selected_cidx = id_to_cidx[closest_match(curr_cluster, hit_id)];

  var selections = d3.selectAll("#svg-timetree > svg > rect:not(.clickedH)")
    .filter(function() {
      return hit_id.includes(parseInt(this.id.substring(3)))
    });

  var cluster = select_cluster(selected_cidx);

  d3.select("#svg-timetree")
    .selectAll("rect:not(.clicked):not(.clickedH)")
    .attr("class","not_SelectedCluster");

  for (const node of selections.nodes()){
    d3.select(node).attr("class", "SelectedCluster");
  }
  cluster.attr("class", "SelectedCluster clicked");
  beadplot(cluster.datum().cluster_idx);

  // Beads in Cluster
  points_ui = d3.selectAll("#svg-cluster > svg > g > circle")
    .filter(function(d) {
      return bead_id_to_accession.includes(d.accessions[0])
    });

  select_beads(points_ui);
  var current_bead_id = bead_hits[cluster_hits[selected_cidx]];
  var total_points = Object.keys(bead_hits).length;

 // Update the search resutls array with the hits
 const stats = search_results.update({
  beads: bead_hits,
  clusters_first_bead: cluster_hits,
  clusters_last_bead: cluster_hits_last_id,
  current_point: current_bead_id,
  total_points: total_points,
  hit_ids: hit_id,
 });

 update_search_stats(stats);

}

/**
 * This function searches for the closest cluster with a search hit (Binary search)
 * @param non_hit_cluster_index: cidx value
 * @param hit_ids: Sorted list of ids (not cidx) of clusters with a hit for the search query
 */
function closest_match(non_hit_cluster_index, hit_ids) {
  const target_index = map_cidx_to_id[non_hit_cluster_index];

  if (target_index <= hit_ids[0])
    return hit_ids[0];
  
  if (target_index >= hit_ids[hit_ids.length - 1])
    return hit_ids[hit_ids.length - 1];

  var i = 0, j = hit_ids.length, mid = 0;
  while (i < j) {
    mid = Math.floor((i + j)/2);
    if (hit_ids[mid] == target_index)
      return hit_ids[mid];

    if (target_index < hit_ids[mid]) {
      if (mid > 0 && target_index > hit_ids[mid - 1])
        return getClosest(hit_ids[mid-1], hit_ids[mid], target_index);
      j = mid;
    }
    else {
      if (mid < hit_ids.length - 1 && target_index < hit_ids[mid+1])
        return getClosest(hit_ids[mid], hit_ids[mid+1], target_index)
      i = mid + 1;
    }
  }
}

function getClosest(key1, key2, target) {
  if (target - key1 >= key2 - target)
    return key2;
  else
    return key1;
}

/**
 * This function is similar to the closest_match function but returns the greater id value (previous cluster)
 * @param {String} non_hit_cluster_index 
 * @param {Array} hit_ids 
 */
function previous_closest_match(non_hit_cluster_index, hit_ids) {
  const target_index = map_cidx_to_id[non_hit_cluster_index];

  if (target_index <= hit_ids[0])
    return hit_ids[0];
  
  if (target_index >= hit_ids[hit_ids.length - 1])
    return hit_ids[hit_ids.length - 1];

  var i = 0, j = hit_ids.length, mid = 0;
  while (i < j) {
    mid = Math.floor((i + j)/2);
    if (hit_ids[mid] == target_index)
      return hit_ids[mid];

    if (target_index < hit_ids[mid]) {
      if (mid > 0 && target_index > hit_ids[mid - 1])
        return hit_ids[mid];
      j = mid;
    }
    else {
      if (mid < hit_ids.length - 1 && target_index < hit_ids[mid+1])
        return hit_ids[mid+1];
      i = mid + 1;
    }
  }
}

/**
 * This function handles search by lineage
 * @param {String} text_query 
 */
function lineage_search(text_query, start_date, end_date) {
  var cidx = lineage_to_cid[text_query.toUpperCase()];

  // Terminates if there is no match
  if (cidx === undefined) {
    $('#error_message').text(`No matches. Please try again.`);
    return;
  }

  var cluster = select_cluster("cidx-"+cidx);

  search_hits = beaddata[cidx].points.filter(function(bead) {
	  temp = (bead.x >= start_date && bead.x <= end_date);
	  return temp;
  });

  d3.select("#svg-timetree")
    .selectAll("rect:not(.clicked):not(.clickedH)")
    .attr("class","not_SelectedCluster");

  cluster.attr("class", "SelectedCluster clicked");
  beadplot(cluster.datum().cluster_idx);

  if (search_hits.length === 0) {
    // Generates table for all points
    var cluster_regions = cluster.datum();
    gentable(cluster_regions);
    draw_region_distribution(cluster_regions.allregions);
    gen_details_table(beaddata[cidx].points);  // update details table with all samples
    $('#error_message').text(`No matches. Please try again.`);
    return;
  }
  else {
    var bead_hits = search_hits.reduce(function(map, bead, i) {
      map[bead.accessions[0]] = i; // value for the search indexing
      return map;
    }, {});

    var hit_id = [map_cidx_to_id["cidx-"+cidx]];
    var cluster_hits = [], cluster_hits_last_id = [];
    var bead_data = beaddata[cidx].points;
    cluster_hits["cidx-"+cidx] = bead_data[0].accessions[0];
    cluster_hits_last_id["cidx-"+cidx] = bead_data[bead_data.length - 1].accessions[0];
    var working_bead = d3.selectAll('circle[id="'+cluster_hits["cidx-"+cidx]+'"]').nodes()[0];
    working_bead.scrollIntoView({block: "center"});
    update_table_individual_bead(d3.select(working_bead).datum());

    const stats = search_results.update({
      beads: bead_hits,
      clusters_first_bead: cluster_hits,
      clusters_last_bead: cluster_hits_last_id,
      current_point: 0,
      total_points: bead_data.length,
      hit_ids: hit_id,
     });
    
     update_search_stats(stats);
  }
  
  // const stats = search_results.update({
  //   current_point: 0,
  //   total_points: 1,
  //  });
  
  // update_search_stats(stats); 
}

/**
 * This function handles search by accession
 * @param {String} text_query 
 */
function accession_search(text_query) {
  var cidx = accn_to_cid[text_query.toUpperCase()];

  if (cidx === undefined) {
    $('#error_message').text(`No matches. Please try again.`);
    return;
  }

  var cluster = select_cluster("cidx-"+cidx);

  d3.select("#svg-timetree")
    .selectAll("rect:not(.clicked):not(.clickedH)")
    .attr("class","not_SelectedCluster");

  cluster.attr("class", "SelectedCluster clicked");
  
  beadplot(cluster.datum().cluster_idx);

  const stats = search_results.update({
    current_point: 0,
    total_points: 1,
   });
  
  update_search_stats(stats); 

  d3.selectAll("circle:not(.selectionH)")
    .attr("class", "not_SelectedBead");

  d3.select("div#svg-cluster")
    .selectAll("line")
    .attr("stroke-opacity", 0.3);

  var bead = d3.selectAll("circle").filter(function(d) {
    return d.accessions.includes(text_query.toUpperCase());
  });

  create_selection(bead);
  bead.node().scrollIntoView({block: "center"});

  update_table_individual_bead(bead.datum());
}

/**
 * This function wraps the main search function. 
 * It reads the search arguments from the ui and ensures that the arguments required by main_search are corectly formatted.
 * It also uses the results produced by the main_search function to populate the ui with the search results
 * When the user changes the search query this should be rerun.
 */
function wrap_search() {
  
  var start_date_text = $('#start-date').val();
  var end_date_text = $('#end-date').val();
  var query = $('#search-input').val();

  var start_date = (start_date_text === "") ? new Date("2020-01-01") : new Date(start_date_text);
  var end_date = (end_date_text === "") ? new Date() : new Date(end_date_text);

  // Exception handing
  if (start_date > end_date) {
    $('#error_message').text(`Start Date must be before the End Date.`);
    return;
  }

  if (isAccn(query)) 
    accession_search(query);
  else if (isLineage(query))
    lineage_search(query, start_date, end_date);
  else 
    main_search(beaddata, query, start_date, end_date);
}


function select_cluster(cidx) {
  d3.selectAll("rect.clicked").attr('class', "default");
  d3.selectAll("rect.clickedH").remove();

  d3.selectAll("circle").dispatch("mouseout");
  var cluster = d3.selectAll('rect[cidx="'+cidx+'"]');

  draw_cluster_box(cluster);
  cluster.nodes()[0].scrollIntoView({block: "center"});

  return cluster;
}

function update_table_individual_bead(bead) {
  draw_halo(bead);
  gentable(bead);
  draw_region_distribution(tabulate(bead.region));
  gen_details_table(bead);
}

function select_next_prev_bead(bead_id_to_accession, curr_bead) {
  d3.selectAll('rect[class="clicked"]').attr('class', "not_SelectedCluster");

  var next_cluster = d3.selectAll('rect[cidx="cidx-'+accn_to_cid[bead_id_to_accession[curr_bead]]+'"]');
  d3.selectAll("rect.clickedH").remove();
  d3.selectAll(".SelectedCluster.clicked").attr('class', 'SelectedCluster'); 
  next_cluster.attr("class", "SelectedCluster clicked");
  draw_cluster_box(next_cluster);
  next_cluster.nodes()[0].scrollIntoView({block: "center"});

  beadplot(next_cluster.datum().cluster_idx);

  // Beads in Cluster
  points_ui = d3.selectAll("#svg-cluster > svg > g > circle")
  .filter(function(d) {
    return bead_id_to_accession.includes(d.accessions[0])
  });
  selected = points_ui.nodes();

  d3.selectAll("circle:not(.selectionH)")
      .attr("class", "not_SelectedBead");

  d3.select("div#svg-cluster")
      .selectAll("line")
      .attr("stroke-opacity", 0.3);
  for (const node of selected) {
    selected_obj = d3.select(node);
    create_selection(selected_obj);
  }
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
    update_table_individual_bead(d3.select(selected[0]).datum())
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
			.map(x => x[1]);
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

function index_lineage(clusters) {
  var index = {};
  for (const cid in clusters) {
    var accns = clusters[cid].lineage
    index[accns] = cid;
  }
  return index;
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
