function index_accession(clusters){
	var keys = [];
	for(cid in clusters) {
		var childkeys =[];
		for(var k in clusters[cid].nodes){ 
			childkeys.push(k) 
			
		}
		keys.push(childkeys);
	}
	return(keys);
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