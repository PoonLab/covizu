// to store references to SVG objects (nodes)
var selected = [];

function select_beads_by_substring(substr) {
	selected = [];
	// this only works for the currently displayed SVG
	d3.selectAll("#svg-cluster > svg > g > circle").filter(function(d) {
		return(d.labels.some(x => x.includes(substr)));
	});
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
	const prefix = /^(E|I|EP|IS|EPI_?|ISL_?)$/i;
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

function search() {
	select_bead_by_accession($('#search-input').val());
}

$('#search-button').on('click', search);
