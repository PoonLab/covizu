// to store references to SVG objects (nodes)
var selected = [];

/**
 * Highlight clusters containing samples that match substring,
 * and highlight samples in currently displayed bead plot.
 * When user switches bead plots, update highlighted samples.
 *
 * @param substr
 */
function select_beads_by_substring(substr) {
	selected = [];
	d3.selectAll("circle").dispatch('mouseout');

	if (substr === "") {
		// user submitted empty string
		clear_selection();
		return;
	}

	var rects = d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)").filter(function(d) {
		return(d.searchtext.match(substr) !== null);
	});

	// jump to the first hit if new search
	if (d3.selectAll('.clicked').empty()){
	  
	  var first_cluster_idx;
	  d3.select(rects.nodes().pop()).attr("class", "clicked");
	  
	  d3.select(rects.nodes().pop()).each(function(d) {
	  first_cluster_idx = d.cluster_idx;
	  create_clusterH(d, this);

	  });
	  
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
	
	var beads = d3.selectAll("#svg-cluster > svg > g > circle").filter(function(d) {
		return(d.labels.some(x => x.includes(substr)));
	});
	selected = beads.nodes();
	d3.selectAll("circle:not(.selectionH)").attr("class", "not_SelectedBead");
	d3.select("div#svg-cluster").selectAll("line").attr("stroke-opacity", 0.3);
	for (const node of beads.nodes()) {
		//d3.select(node).dispatch('mouseover');
		var selected_obj = d3.select(node);
		create_selection(selected_obj);
	}
	if (beads.nodes().length !== 0) {
	  beads.nodes()[0].scrollIntoView({block: "center"});
	}
}


/**
 * Highlight and jump to node in beadplot by sample accession number.
 * @param {string} accn:  accession number to search for
 */
function select_bead_by_accession(accn) {
  
  d3.selectAll("circle").dispatch('mouseout');
  // switch to cluster beadplot
  var cid = accn_to_cid[accn];
  
  if (cid !== undefined) {
     if (d3.selectAll('.clicked').empty()) {
       
        d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)").attr("class", "not_SelectedCluster");
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
        
     } else {
       d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)").attr("class", "not_SelectedCluster");
       var rect = d3.selectAll("#svg-timetree > svg > g > rect:not(.clickedH)")
        .filter(function(d) { return(d.cluster_idx === cid); })
        .attr("class", "SelectedCluster");
        
        var bead = d3.selectAll("circle").filter(function(d) {
          return d.accessions.includes(accn);
        });
        
        if (bead.nodes().length !== 0) {
          create_selection(bead);
          bead.node().scrollIntoView({block: "center"});
        } else {
          d3.select("div#svg-cluster").selectAll("line").attr("stroke-opacity", 0.3);
          
          d3.select("div#svg-cluster")
          .selectAll("circle:not(.SelectedBead):not(.selectionH)")
          .attr("class", "not_SelectedBead");
        }
     }
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

function search() {
	var query = $('#search-input').val();
	const accn_pat = /^EPI_ISL_[0-9]+$/i;  // case-insensitive
	if (accn_pat.test(query)) {
		// user is querying accession number
		select_bead_by_accession(query);
	}
	else {
		// substring search
		select_beads_by_substring(query);
	}
}

$('#search-input').on('keydown', function(e) {
	if (e.keyCode == 13) {
		// type <enter> to run search
		if ($('#search-input').val() !== "") {
		  d3.selectAll("rect.clicked").attr('class', "default");
		}
		d3.selectAll("rect.clickedH").remove();
		search();
	}
})