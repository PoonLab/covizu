const compression = require('compression');
const express = require('express');
const app = express();
const clusters = require('./data/clusters.json')
const { utcDate } = require('./server/utils')
const { parse_clusters, map_clusters_to_tips, index_accessions, index_lineage, get_recombinants } = require('./server/parseCluster')
const { readTree } = require('./server/phylo')
const fs = require('fs');
require('dotenv').config();

var http = require('http');
var https = require('https');

if (process.env.PROD) {
  var credentials = {
    key: fs.readFileSync(process.env.PRVTKEY), 
    cert: fs.readFileSync(process.env.CRT)
  };
}

app.use(compression());

try {
  var tree = fs.readFileSync('./data/timetree.nwk', 'utf8');
} catch(e) {
  console.log('Error:', e.stack);
}

const df = readTree(tree)
const beaddata = parse_clusters(clusters)
const { tips, recombinant_tips } = map_clusters_to_tips(df, clusters)
const accn_to_cid = index_accessions(clusters)
const lineage_to_cid = index_lineage(clusters)
const recombinants = get_recombinants()

// This is a hack to match anything that could be an acc number prefix
const prefix = /^(E|I|EP|IS|EPI_I|EPI_IS|EPI_ISL_?|EPI_?|ISL_?|[A-Z]\.[1-9]+)$/i;
const MIN_RESULTS = 10;
const normalize = (str) => str.replace(/[^a-z0-9]/gi, '').toLowerCase();
const data = Object.keys(accn_to_cid).sort().concat(Object.keys(lineage_to_cid).sort()).map(accn => [
  normalize(accn), accn
]);

app.get('/api/edgeList/:cindex', (req, res) => {
  res.send(beaddata[req.params.cindex].edgelist)
});

app.get('/api/points/:cindex', (req, res) => {
  res.send(beaddata[req.params.cindex].points)
});

app.get('/api/variants/:cindex', (req, res) => {
  res.send(beaddata[req.params.cindex].variants)
});

app.get('/api/tips', (req, res) => {
  res.send(tips)
});

app.get('/api/df', (req, res) => {
  res.send(df)
});

app.get('/api/lineage/:cindex', (req, res) => {
  res.send(clusters[req.params.cindex].lineage)
});

app.get('/api/cid/:accession', (req, res) => {
  res.send(accn_to_cid[req.params.accession])
});

app.get('/api/cid', (req, res) => {
  res.send(accn_to_cid)
});

app.get('/api/lineagetocid', (req, res) => {
  res.send(lineage_to_cid)
});

app.get('/api/recombtips', (req, res) => {
  res.send(recombinant_tips)
});

app.get('/api/searchHits/:query/:start/:end', (req, res) => {
  // Flatten the json data to an array with bead data only
  let flat_data = beaddata.map(bead => bead.points).flat();
  let start_date = utcDate(req.params.start);
  let end_date = utcDate(req.params.end)

  //Find all the beads that are a hit. Convert text_query to lower case and checks to see if there is a match
  let search_hits = flat_data.filter(function(bead) {
	  temp = (bead.accessions.some(accession => (accession.toLowerCase()).includes(req.params.query.toLowerCase())) || 
		  bead.labels.some(label => (label.toLowerCase()).includes(req.params.query.toLowerCase()))) && 
		  (bead.x >= start_date && bead.x <= end_date);
	  return temp;
  });

  res.send(search_hits)
});

app.get('/api/getHits/:query', (req, res) => {
  function as_label(search_data) {
    const [, accn] = search_data;
    return accn;
  }

  const term = req.params.query

  if (!/\d/.test(term)) {
    if (prefix.test(term)) {
      res.send(data.slice(0, MIN_RESULTS).map(as_label));
    } else {
      res.send([]);
    }
  } else {
    const result = data.filter(array => array[0].indexOf(normalize(term)) > -1);
    res.send(result.slice(0, MIN_RESULTS).map(as_label));
  }
});

const port = process.env.PORT || 8001;

app.use(express.static('.'));

// For the prod environment, need to create a https server
if (process.env.PROD) {
  var httpsServer = https.createServer(credentials, app);
  httpsServer.listen(8002, () => console.log(`Listening on Port ${8002}...`))
}

var httpServer = http.createServer(app);
httpServer.listen(port, () => console.log(`Listening on Port ${port}...`));
