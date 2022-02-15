const express = require('express');
const app = express();
const clusters = require('./data/clusters.json')
const { utcDate } = require('./utils/helpers')
const { parse_clusters, readTree, map_clusters_to_tips, index_accessions, index_lineage } = require('./utils/parseCluster')

const fs = require('fs');

try {
    var tree = fs.readFileSync('./data/timetree.nwk', 'utf8');
    // console.log(tree);    
} catch(e) {
    console.log('Error:', e.stack);
}

const df = readTree(tree)
const beaddata = parse_clusters(clusters)
const tips = map_clusters_to_tips(df, clusters)
const accn_to_cid = index_accessions(clusters)
const lineage_to_cid = index_lineage(clusters)

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
})

app.get('/api/cid/:accession', (req, res) => {
  res.send(accn_to_cid[req.params.accession])
})

app.get('/api/cid', (req, res) => {
  res.send(accn_to_cid)
});

app.get('/api/lineagetocid', (req, res) => {
  res.send(lineage_to_cid)
})

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
})

const port = process.env.PORT || 8001;

app.listen(port, () => console.log(`Listening on Port ${port}...`))
app.use(express.static('.'));