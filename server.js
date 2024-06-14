const { $HTTP_PORT, $NODE_ENV } = require("./config/config")
const {
  $COVIZU_CONNECTION_URI, $ACTIVE_DATABASE, $COLLECTION__CLUSTERS
} = require('./config/dbconfig')
require("./globalVariables");

const compression = require('compression');
const express = require('express');
const app = express();
const { utcDate } = require('./server/utils')

var http = require('http');

app.use(compression());
app.use(express.static('.'));

let dbManager;
if ($NODE_ENV == 'TEST') {
  const { TestDBManager } = require("./dbmanager_test");
  dbManager = new TestDBManager();
  dbManager.load_globalData();
}
if ($NODE_ENV == 'DEV' || $NODE_ENV == 'PROD') {
  const { DBManager } = require("./dbmanager");
  dbManager = new DBManager();
  dbManager.set_dbUrl($COVIZU_CONNECTION_URI)
  dbManager.set_dbName($ACTIVE_DATABASE)
  dbManager.start()
}

app.get('/api/edgeList/:cindex', (req, res) => {
  dbManager.get_edgeList(req.params.cindex).then(result=>{
    res.send(result);
  })
});

app.get('/api/points/:cindex', (req, res) => {
  dbManager.get_points(req.params.cindex).then(result=>{
    res.send(result);
  })
});

app.get('/api/variants/:cindex', (req, res) => {
  dbManager.get_variants(req.params.cindex).then(result=>{
    res.send(result);
  })
});

app.get('/api/lineage/:cindex', (req, res) => {
  dbManager.get_lineage(req.params.cindex).then(result=>{
    res.send(result);
  })
});

app.get('/api/tips', (req, res) => {
  dbManager.get_tips().then(result=>{
    res.send(result);
  })
});

app.get('/api/display/:cindex', (req, res) => {
  dbManager.get_lineage(req.params.cindex).then(result => {
    dbManager.get_display(result).then(display => {
      res.send(display);
    })
  })
});

app.get('/api/df', (req, res) => {
  dbManager.get_df().then(result=>{
    res.send(result)
  })
});

app.get('/api/xbb', (req, res) => {
  dbManager.get_xbb_df().then(result=>{
    res.send(result)
  })
});

app.get('/api/regionmap', (req, res) => {
  dbManager.get_regionMap().then(result=>{
    res.send(result);
  })
  // res.send(dbManager.get_regionMap());
})

app.get('/api/cid/:accession', (req, res) => {
  dbManager.get_accession(req.params.accession.toUpperCase()).then(result=>{
    res.send(result);
  })
});

// this endpoint serves a massive file ~60MB. should we keep it?
app.get('/api/cid', (req, res) => {
  dbManager.get_cid().then(result=>{
    res.send(result);
  })
  // res.send(dbManager.get_cid());
});

app.get('/api/lineagetocid', (req, res) => {
  dbManager.get_lineageToCid().then(result=>{
    res.send(result);
  })
  // res.send(dbManager.get_lineageToCid());
});

app.get('/api/recombtips', (req, res) => {
  dbManager.get_recombinantTips().then(result=>{
    res.send(result);
  })
  // res.send(dbManager.get_recombinantTips());
});

app.get('/api/searchHits/:query/:start/:end', (req, res) => {
  const start_date = utcDate(req.params.start);
  const end_date = utcDate(req.params.end);
  const query = req.params.query.toLocaleLowerCase();
  dbManager.get_searchHits(start_date, end_date, query).then(result=>{
    res.send(result);
  })
  // res.send(dbManager.get_searchHits(start_date, end_date, query));
});

app.get('/api/getHits/:query', (req, res) => {
  const term = req.params.query
  dbManager.get_hits(term).then(result=>{
    res.send(result);
  })
  // res.send(dbManager.get_hits(term));
});

var httpServer = http.createServer(app);
httpServer.listen($HTTP_PORT, () => console.log(`Listening on Port ${$HTTP_PORT}...`));
