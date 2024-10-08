const { MongoClient } = require("mongodb");
const readline = require("readline");
const fs = require('fs');
const { $NODE_ENV } = require("../config/config")

const rl = readline.createInterface({
  input: process.stdin,
  output: process.stdout
});

const {
  parse_clusters,
  map_clusters_to_tips,
  index_accessions,
  index_lineage
} = require('../server/parseCluster')

const { readTree } = require("../server/phylo");

// we have to use globalVariables because the parse_clusters updates the global.region_map on the fly.
require("../globalVariables")

const {
  $ADMIN_CONNECTION_URI,
  $ACTIVE_DATABASE,
  $COLLECTION__CLUSTERS,
  $COLLECTION__DBSTATS,
  $COLLECTION__BEADDATA,
  $COLLECTION__TIPS,
  $COLLECTION__RECOMBINANT_TIPS,
  $COLLECTION__ACCN_TO_CID,
  $COLLECTION__LINEAGE_TO_CID,
  $COLLECTION__REGION_MAP,
  $COLLECTION__DF_TREE,
  $COLLECTION__XBB_TREE,
  $COLLECTION__AUTOCOMPLETE_DATA,
  $COLLECTION__FLAT_DATA,
  $PROJECT_ROOT,
  $JSON_DATA_FOLDER,
  $JSONFILE__CLUSTERS,
  $JSONFILE__DBSTATS,
  $NWKFILE__TREE,
  $XBB__TREE
} = require('../config/dbconfig')


const adminUserConnection = new MongoClient($ADMIN_CONNECTION_URI);
console.log(
  `ADMIN_CONNECTION_URI:            ${$ADMIN_CONNECTION_URI}
   ACTIVE_DATABASE:                 ${$ACTIVE_DATABASE}
   COLLECTION__CLUSTERS:            ${$COLLECTION__CLUSTERS}
   COLLECTION__BEADDATA:            ${$COLLECTION__BEADDATA}
   COLLECTION__TIPS:                ${$COLLECTION__TIPS}
   COLLECTION__RECOMBINANT_TIPS     ${$COLLECTION__RECOMBINANT_TIPS}
   COLLECTION__ACCN_TO_CID          ${$COLLECTION__ACCN_TO_CID}
   COLLECTION__LINEAGE_TO_CID       ${$COLLECTION__LINEAGE_TO_CID}
   COLLECTION__REGION_MAP           ${$COLLECTION__REGION_MAP}
   COLLECTION__DF_TREE              ${$COLLECTION__DF_TREE}
   COLLECTION__AUTOCOMPLETE_DATA    ${$COLLECTION__AUTOCOMPLETE_DATA}
   COLLECTION__FLAT_DATA            ${$COLLECTION__FLAT_DATA}
   PROJECT_ROOT                     ${$PROJECT_ROOT}
   JSON_DATA_FOLDER                 ${$JSON_DATA_FOLDER}
   JSONFILE__CLUSTERS               ${$JSONFILE__CLUSTERS}
   JSONFILE__DBSTATS                ${$JSONFILE__DBSTATS}
   NWKFILE_TREE                     ${$NWKFILE__TREE}
   XBB__TREE                        ${$XBB__TREE}`
)

if(!$ACTIVE_DATABASE)
{
    throw new Error("Environment variable DBNUMBER must be 1 or 2");
}

async function updateDatabase() {

  console.log('starting connection...');
  adminUserConnection.connect()
    .then((res) => {
      console.log("Connecting as covizu_admin...")
    })
    .catch((error) => {
      console.log(`Error while connecting as covizu_admin:${error}`);
      throw new Error(error);
    })
    .then(() => {
      console.log("CONNECTED!");
    })
    .then(async () => {
      let db;

      db = adminUserConnection.db($ACTIVE_DATABASE);
      
      /** delete all collections */
      const existingCollections = await db.listCollections().toArray();
      existingCollections.forEach((e)=>{
        console.log(`Deleting collection ${$ACTIVE_DATABASE}.${e.name}`);
        db.collection(e.name).drop();
      })
      
      let textfile;
      textfile = $PROJECT_ROOT+"/"+$JSON_DATA_FOLDER+"/"+$NWKFILE__TREE;
      console.log(`Reading ${textfile}`);
      try {
        // global.tree = fs.readFileSync(`${$JSON_DATA_FOLDER}/${$NWKFILE__TREE}`, 'utf8');
        global.tree = fs.readFileSync(textfile, 'utf8');
      }
      catch (e) {
        console.error(`Failed reading ${textfile} : `, e);
        return;
      }

      textfile = $PROJECT_ROOT+"/"+$JSON_DATA_FOLDER+"/"+$XBB__TREE;
      console.log(`Reading ${textfile}`);
      try {
        global.xbbtree = fs.readFileSync(textfile, 'utf8');
      }
      catch (e) {
        console.error(`Failed reading ${textfile} : `, e);
        return;
      }

      textfile = $PROJECT_ROOT+"/"+$JSON_DATA_FOLDER+"/"+$JSONFILE__CLUSTERS;
      console.log(`Reading ${textfile}`);
      try{
        global.clusters = require(textfile);
      }
      catch(e)
      {
        console.error(`Failed reading ${textfile} : `, e);
      }

      textfile = $PROJECT_ROOT+"/"+$JSON_DATA_FOLDER+"/"+$JSONFILE__DBSTATS;
      console.log(`Reading ${textfile}`);
      try{
        global.dbstats = require(textfile);
      }
      catch(e)
      {
        console.error(`Failed reading ${textfile} : `, e);
      }
      
      console.log("Preparing beaddata from clusters");
      global.beaddata = parse_clusters(global.clusters);

      console.log("Preparing df from tree");
      global.df = readTree(global.tree);
      global.df_xbb = readTree(global.xbbtree);

      console.log("Preparing tips, recombinant_tips from df,clusters")
      const { tips, tips_xbb, recombinant_tips } = map_clusters_to_tips(global.df, global.df_xbb, global.clusters);
      global.tips = tips;
      global.recombinant_tips = recombinant_tips;

      console.log("Preparing accn_to_cid from clusters")
      global.accn_to_cid = index_accessions(global.clusters)

      console.log("Preparing lineage_to_cid from clusters");
      global.lineage_to_cid = index_lineage(global.clusters)

      console.log("Preparing autocomplete_data from accn_to_cid and lineage_to_cid");
      global.autocomplete_data = Object.keys(global.accn_to_cid).sort()
      .concat(Object.keys(global.lineage_to_cid).sort())
      .map(accn => [global.normalize(accn), accn]);

      console.log("Preparing flat_data from beaddata");
      global.flat_data = global.beaddata.map(bead=>bead.points).flat();
      
      
      console.log("Writing to database...")
      let res;

      res = await db.collection($COLLECTION__CLUSTERS).insertMany(global.clusters);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__CLUSTERS}`);
      delete global.clusters;

      global.lineage_stats = Object.entries(global.dbstats['lineages']).map(([key, value]) => ({ _id: key, ...value}));
      res = await db.collection($COLLECTION__DBSTATS).insertMany(global.lineage_stats);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__DBSTATS}`);
      delete global.dbstats;
      delete global.lineage_stats;

      res = await db.collection($COLLECTION__BEADDATA).insertMany(global.beaddata);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__BEADDATA}`);
      delete global.beaddata;

      res = await db.collection($COLLECTION__TIPS).insertMany(global.tips);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__TIPS}`);
      delete global.tips;

      if (global.recombinant_tips.length > 0) {
        res = await db.collection($COLLECTION__RECOMBINANT_TIPS).insertMany(global.recombinant_tips);
        console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__RECOMBINANT_TIPS}`);
        delete global.recombinant_tips;
      }
      
      /** MongoDB has a size-limit of ~17MB per record. So accn_to_cid needs to be broken down into an array of objects*/
      res = await db.collection($COLLECTION__ACCN_TO_CID).insertMany(Object.entries(global.accn_to_cid).map(el => { j = {}; j[el[0]] = el[1]; return j }));
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__ACCN_TO_CID}`);
      delete global.accn_to_cid;

      /** MongoDB has a size-limit of ~17MB per record. So lineage_to_cid needs to be broken down into an array of objects */
      res = await db.collection($COLLECTION__LINEAGE_TO_CID).insertMany(Object.entries(global.lineage_to_cid).map(el => { j = {}; j[el[0]] = el[1]; return j }));
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__LINEAGE_TO_CID}`);
      delete global.lineage_to_cid;

      res = await db.collection($COLLECTION__REGION_MAP).insertOne(global.region_map);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__REGION_MAP}`);
      delete global.region_map;

      res = await db.collection($COLLECTION__DF_TREE).insertMany(global.df);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__DF_TREE}`);
      delete global.df;

      res = await db.collection($COLLECTION__XBB_TREE).insertMany(global.df_xbb);
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__DF_TREE}`);
      delete global.df_xbb;

      res = await db.collection($COLLECTION__AUTOCOMPLETE_DATA).insertMany(global.autocomplete_data.map(subArray => {return {'norm': subArray[0],'accn': subArray[1]};}))
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__AUTOCOMPLETE_DATA}`);
      delete global.autocomplete_data;

      res = await db.collection($COLLECTION__FLAT_DATA).insertMany(global.flat_data)
      console.log(`Created ${res.insertedCount} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__FLAT_DATA}`);
      delete global.flat_data;


      res = await adminUserConnection.close();
      console.log("Sucessfully wrote to database and closed connection!");

    })
    .catch(error => {
      console.log(`ERROR: ${error}`)
      adminUserConnection.close().then(() => {
        console.log("Failed to write to database. Closing connection.");
      })
    })
}

const question_warning = "This will delete records within database " + $ACTIVE_DATABASE + ". Are you sure? (Y/n)";
const userConfirmationCallback = function (answer) {
  if (answer == "Y") {
    rl.close();
    updateDatabase();
  }
  else if (answer == 'n') {
    rl.close();
    console.log("EXITING");
    process.exit();
  }
  else {
    rl.question(question_warning, userConfirmationCallback);
  }
}

if ($NODE_ENV == 'PROD' || $NODE_ENV == 'DEV') {
  rl.close()
  updateDatabase();
}
else {
  rl.question(question_warning, userConfirmationCallback)
}