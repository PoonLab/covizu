// to store references to SVG objects (nodes)
var selected = [];

/**
* Prepare the search stats with initial stats
*/
function prepare_search_stats(initial_stats) {
  let stats = {
    ...initial_stats,
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
  start_idx: [],
});

function update_search_stats(stats) {
  $('#search_stats').text(`${stats.current_point} of ${stats.total_points} points`);
}

function find_beads_points(beadsdata){
  return beadsdata.map(bead => bead.points).flat()
}

/**
 * Highlight clusters for which the cluster function returns true
 * and highlight points in currently displayed bead plot.
 * When user switches bead plots, update highlighted points.
 *
 * @param clusters function that takes a cluster datum and returns true
 *                 when the cluster must be selected. False otherwise.
 * @param points function that takes a points datum and returns true
 *               when the point must be selected. False otherwise.
 */
function select_beads_by_comparators(clusters, points) {
  selected = [];
	d3.selectAll("circle").dispatch('mouseout');

	var rects = d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)").filter(clusters);

	// jump to the first hit if new search
	if (d3.selectAll('.clicked').empty()){

	  var first_cluster_idx;
	  d3.select(rects.nodes().pop()).attr("class", "clicked");
	  d3.select(rects.nodes().pop()).each(function(d) {
      first_cluster_idx = d.cluster_idx;
  	  create_clusterH(d, this);
	  });
    first_cluster_idx = typeof first_cluster_idx == "undefined" ? Math.floor(Math.random() * 239) : first_cluster_idx
    beadplot(first_cluster_idx);

	  var itter = rects.nodes().splice(0, rects.nodes().length-1);

	  d3.select("#svg-timetree").selectAll("rect:not(.clicked):not(.clickedH)").attr("class","not_SelectedCluster");

	} else{

	  d3.select("#svg-timetree").selectAll("rect:not(.clickedH)").attr("class","not_SelectedCluster");

	  var itter = rects.nodes();
	}

	// FIXME: other matching clusters not getting highlighted, due to click below?
	for (const node of itter) {
		d3.select(node).attr("class","SelectedCluster");
	}

	var points_ui = d3.selectAll("#svg-cluster > svg > g > circle").filter(points);
	selected = points_ui.nodes();
	d3.selectAll("circle:not(.selectionH)").attr("class", "not_SelectedBead");
	d3.select("div#svg-cluster").selectAll("line").attr("stroke-opacity", 0.3);
	for (const node of points_ui.nodes()) {
		//d3.select(node).dispatch('mouseover');
		var selected_obj = d3.select(node);
		create_selection(selected_obj);
	}
	if (points_ui.nodes().length !== 0) {
	  points_ui.nodes()[0].scrollIntoView({block: "center"});
	}
}

/**
 * Highlight clusters containing samples that match substring,
 * and highlight samples in currently displayed bead plot.
 * When user switches bead plots, update highlighted samples.
 *
 * @param substr
 */
function select_beads_by_substring(substr) {
  if (substr === "") {
    // user submitted empty string
    clear_selection();
    return;
  }

  const clusters = (datum) => datum.searchtext.match(substr) !== null;
  const points = (datum) => datum.labels.some(x => x.includes(substr));
  select_beads_by_comparators(clusters, points);
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
      create_clusterH(d, this);
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
    draw_region_distribution(table(datum.region));
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

function search(beaddata) {
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

function search_by_dates(beaddata, start, end){
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
}
