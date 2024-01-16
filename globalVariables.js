global.clusters; // mock
global.beaddata; // mock

global.tree;
global.df; // mock
global.tips; // mock
global.recombinant_tips; // mock
global.accn_to_cid; // mock
global.lineage_to_cid; // mock
global.region_map = {}; // mock

global.prefix = /^(E|I|EP|IS|EPI_I|EPI_IS|EPI_ISL_?|EPI_?|ISL_?|[A-Z]\.[1-9]+)$/i;
global.MIN_RESULTS = 10;
global.normalize = (str) => str.replace(/[^a-z0-9]/gi, '').toLowerCase();
global.autocomplete_data; 
global.flat_data; 
