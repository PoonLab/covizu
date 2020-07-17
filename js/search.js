// to store references to SVG objects (nodes)
var selected = [];

/**
 * Highlight clusters containing samples that match substring,
 * and highlight samples in currently displayed bead plot.
 * When user switches bead plots, update highlighted samples.
 *
 * @param substr
 */
function select_beads_by_substring(substr, accn) {
	selected = [];
	d3.selectAll("circle").dispatch('mouseout');

	if (substr === "") {
		// user submitted empty string
		return;
	}

	var rects = d3.selectAll("#svg-timetree > svg > g > rect").filter(function(d) {
		return(d.searchtext.match(substr) !== null);
	});

	// jump to the first hit
	//d3.select(rects.nodes()[0]).dispatch('click');
	
	// FIXME: other matching clusters not getting highlighted, due to click below?
	d3.selectAll("rect").classed("SelectedCluster", false);
	for (const node of rects.nodes()) {
		d3.select(node).classed("SelectedCluster", true);
	}
	d3.select("#svg-timetree").selectAll("rect:not(.SelectedCluster)")
	.attr("fill-opacity", 0.1);
	d3.select("#svg-timetree").selectAll("lines")
	.attr("stroke-opacity", 0.3);
	
	var beads = d3.selectAll("#svg-cluster > svg > g > circle").filter(function(d) {
		return(d.labels.some(x => x.includes(substr)));
	});
	selected = beads.nodes();
	d3.selectAll("circle:not(.selectionH)").attr("class", "not_SelectedBead");
	//beads.nodes()[0].scrollIntoView();
	for (const node of beads.nodes()) {
		//d3.select(node).dispatch('mouseover');
		var selected_obj = d3.select(node);
		create_selection(selected_obj);
	}
}


/**
 * Highlight and jump to node in beadplot by sample accession number.
 * @param {string} accn:  accession number to search for
 */
function select_bead_by_accession(accn) {
	// reset all highlights
	selected = [];
	d3.selectAll("circle").dispatch('mouseout');

	// switch to cluster beadplot
	var cid = accn_to_cid[accn];
	if (cid !== undefined) {
		var rect = d3.selectAll("#svg-timetree > svg > g > rect")
				.filter(function(d) { return(d.cluster_idx === cid); });
		d3.select(rect.node()).dispatch("click");

		var bead = d3.selectAll("circle").filter(function(d) {
			return d.accessions.includes(accn);
		});
		bead.dispatch('mouseover');
		bead.node().scrollIntoView();
		selected.push(bead.node());
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

function search(accn) {
	var query = $('#search-input').val();
	const accn_pat = /^EPI_ISL_[0-9]+$/i;  // case-insensitive
	if (accn_pat.test(query)) {
		// user is querying accession number
		select_bead_by_accession(query);
	}
	else {
		// substring search
		select_beads_by_substring(query, accn);
	}
}

$('#search-button').on('click', search(0));

$('#search-input').on('keydown', function(e) {
	if (e.keyCode == 13) {
		// type <enter> to run search
		search(0);
	}
})