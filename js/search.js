// to store references to SVG objects
var selected = [];


function select_bead_by_accession(accn) {
	// reset all highlights
	selected = [];
	d3.selectAll("circle").dispatch('mouseout');

	// switch to cluster beadplot
	var cid = accn_to_cid[accn];
	if (cid !== undefined) {
		beadplot(cid);
		var bead = d3.selectAll("circle").filter(function(d) {
			return d.accessions.includes(accn);
		});
		bead.dispatch('mouseover');
		selected.push(bead);
	}

}


/**
 * Populate Object with accession-cluster ID as key-value pairs
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

function graphsearch(keys){
	var term = $('#search-input').val();
	for(cid in keys){
		if (keys[cid].includes(term)){
			beadplot(cid);
		}
		}
	} 


$('#search-button').attr("onClick", 'graphsearch(index_accession(clusters));');