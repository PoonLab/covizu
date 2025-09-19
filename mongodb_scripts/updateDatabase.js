const { MongoClient } = require('mongodb');
const readline = require('readline');
const fs = require('fs');
const path = require('path');

const { chain } = require('stream-chain');
const { parser } = require('stream-json');
const { streamArray } = require('stream-json/streamers/StreamArray');
const { streamObject } = require('stream-json/streamers/StreamObject');

const { $NODE_ENV } = require('../config/config');

const {
  parse_clusters,             
  index_accessions,
  index_lineage,
  map_tips
} = require('../server/parseCluster');

const { utcDate } = require('../server/utils');   
const d3 = require('../js/d3');                   
const { readTree } = require('../server/phylo');

require('../globalVariables');

const rl = readline.createInterface({ input: process.stdin, output: process.stdout });

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
} = require('../config/dbconfig');

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
);

if (!$ACTIVE_DATABASE) {
  throw new Error('Environment variable DBNUMBER must be 1 or 2');
}

function pjoin(rel) { return path.join($PROJECT_ROOT, $JSON_DATA_FOLDER, rel); }

function detectTopLevelType(file) {
  const fd = fs.openSync(file, 'r');
  const b = Buffer.alloc(1);
  let pos = 0;
  try {
    while (fs.readSync(fd, b, 0, 1, pos++) === 1) {
      const ch = String.fromCharCode(b[0]);
      if (!/\s/.test(ch)) {
        return ch === '[' ? 'array' : (ch === '{' ? 'object' : 'unknown');
      }
    }
  } finally { try { fs.closeSync(fd); } catch {} }
  throw new Error(`Empty or unreadable JSON file: ${file}`);
}

async function safeInsertMany(col, docs) {
  if (!docs.length) return 0;
  try {
    const res = await col.insertMany(docs, { ordered: false });
    return res.insertedCount ?? Object.keys(res.insertedIds || {}).length ?? docs.length;
  } catch (err) {
    const res = err.result || err.writeErrors?.result;
    if (res && (res.insertedCount || res.insertedIds)) {
      return res.insertedCount ?? Object.keys(res.insertedIds || {}).length ?? 0;
    }
    throw err;
  }
}

function mergeMap(dst, src) {
  if (!src) return;
  for (const k of Object.keys(src)) {
    if (dst[k]) {
      const seen = new Set(dst[k]);
      for (const v of src[k]) if (!seen.has(v)) dst[k].push(v);
    } else {
      dst[k] = Array.isArray(src[k]) ? src[k].slice() : src[k];
    }
  }
}

async function updateDatabase() {
  console.log('starting connection...');
  const client = new MongoClient($ADMIN_CONNECTION_URI);

  await client.connect()
    .then(() => console.log('Connecting as covizu_admin...'))
    .catch((error) => { console.log(`Error while connecting as covizu_admin:${error}`); throw error; });

  console.log('CONNECTED!');
  const db = client.db($ACTIVE_DATABASE);

  // Drop existing collections
  const existingCollections = await db.listCollections().toArray();
  for (const e of existingCollections) {
    console.log(`Deleting collection ${$ACTIVE_DATABASE}.${e.name}`);
    try { await db.collection(e.name).drop(); } catch {}
  }

  // Get trees
  const treePath = pjoin($NWKFILE__TREE);
  console.log(`Reading ${treePath}`);
  try { global.tree = fs.readFileSync(treePath, 'utf8'); }
  catch (e) { console.error(`Failed reading ${treePath} : `, e); return; }

  const xbbPath = pjoin($XBB__TREE);
  console.log(`Reading ${xbbPath}`);
  try { global.xbbtree = fs.readFileSync(xbbPath, 'utf8'); }
  catch (e) { console.error(`Failed reading ${xbbPath} : `, e); return; }

  // Prepare df
  console.log('Preparing df from tree');
  global.df = readTree(global.tree);
  global.df_xbb = readTree(global.xbbtree);

  // Read dbstats
  const dbstatsPath = pjoin($JSONFILE__DBSTATS);
  console.log(`Reading ${dbstatsPath}`);
  try { global.dbstats = require(dbstatsPath); }
  catch (e) { console.error(`Failed reading ${dbstatsPath} : `, e); return; }

  // collections
  const colClusters    = db.collection($COLLECTION__CLUSTERS);
  const colBeaddata    = db.collection($COLLECTION__BEADDATA);
  const colFlat        = db.collection($COLLECTION__FLAT_DATA);

  // Reset indexes
  try {
    await colBeaddata.createIndex({ cidx: 1 });     
    await colFlat.createIndex({ cidx: 1, date: 1 });
  } catch (e) {
    console.warn('Index creation warning:', e.message);
  }

  const clustersFile = pjoin($JSONFILE__CLUSTERS);
  console.log(`Streaming ${clustersFile}`);
  const topType = detectTopLevelType(clustersFile);
  const streamer = topType === 'array' ? streamArray() : streamObject();
  
  // Stream clusters.json
  const pipeline = chain([
    fs.createReadStream(clustersFile),
    parser(),
    streamer
  ]);

  // batch sizes
  const BATCH_CLUSTERS   = 200;
  const BATCH_BEADDATA   = 100;
  const BATCH_FLAT       = 1500;
  const BATCH_MAP_CHUNKS = 2000;

  // batch arrays
  let clustersBatch = [];
  let beaddataBatch = [];
  let flatBatch     = [];

  let totalClusters = 0;
  let totalBeaddata = 0;
  let totalFlat     = 0;

  const accn_to_cid = Object.create(null);
  const lineage_to_cid = Object.create(null);

  let tips = global.df.filter(x => x.children.length === 0);
  const tip_labels = tips.map(x => x.thisLabel);
  const recombinant_tips_partial = []; 

  function flushClusters()    { return safeInsertMany(colClusters, clustersBatch).then(c => { totalClusters += c; clustersBatch = []; }); }
  function flushBeaddata()    { return safeInsertMany(colBeaddata, beaddataBatch).then(c => { totalBeaddata += c; beaddataBatch = []; }); }
  function flushFlat()        { return safeInsertMany(colFlat,     flatBatch).then(c => { totalFlat     += c; flatBatch     = []; }); }

  pipeline.on('data', async ({ key, value: cluster }) => {
    pipeline.pause();
    try {
      clustersBatch.push(cluster);
      if (clustersBatch.length >= BATCH_CLUSTERS) await flushClusters();

      const container = (topType === 'array') ? [cluster] : { [key]: cluster };
      const beadArr = parse_clusters(container); 
      const bead = beadArr[0];

      beaddataBatch.push({ variants: bead.variants, edgelist: bead.edgelist, points: bead.points });
      if (beaddataBatch.length >= BATCH_BEADDATA) await flushBeaddata();

      for (const p of bead.points) flatBatch.push(p);
      if (flatBatch.length >= BATCH_FLAT) await flushFlat();

      mergeMap(accn_to_cid, index_accessions(container));
      mergeMap(lineage_to_cid, index_lineage(container));

      if (cluster['nodes'].length === 1) {
      } else {
        const raw_lineage = global.dbstats['lineages'][cluster['lineage']]['raw_lineage'] || '';
        if (raw_lineage.startsWith('XBB')) {
          let tips_xbb = global.df_xbb.filter(x => x.children.length === 0);
          const tip_labels_xbb = tips_xbb.map(x => x.thisLabel);
          const labels = Object.keys(cluster['nodes']);
          const root = tip_labels_xbb.filter(v => v === cluster['lineage'])[0];
          if (root !== undefined) {
            tips_xbb = map_tips(key, labels, root, tips_xbb, tip_labels_xbb, cluster);
          } else {
            console.log('Failed to match XBB cluster index', key, 'to a tip in the XBB tree');
          }
        } else if (raw_lineage.startsWith('X')) {
          let labels = Object.keys(cluster['nodes']);

          // collect coldates
          let coldates = [];
          for (let i = 0; i < labels.length; i++) {
            const label = labels[i];
            const variant = cluster['nodes'][label];
            coldates = coldates.concat(variant.map(x => x[0]));
          }
          coldates.sort();

          const first_date = utcDate(coldates[0]);
          const last_date  = utcDate(coldates[coldates.length - 1]);

          const date_diffs = coldates.map(x => (utcDate(x) - first_date) / (1000 * 3600 * 24));
          const mean_date = Math.round(date_diffs.reduce((a, b) => a + b, 0) / date_diffs.length);

          const times = coldates.map(x => utcDate(x).getTime());
          const origin = 18231;
          const mean_time = times.reduce((x, y) => x + y) / times.length / 8.64e7 - origin;
          const rate = 0.0655342;
          const exp_diffs = rate * mean_time;

          const tip_stats = global.dbstats['lineages'][cluster['lineage']];

          Date.prototype.addDays = function(days) { const d = new Date(this.valueOf()); d.setDate(d.getDate() + days); return d; };

          recombinant_tips_partial.push({
            cluster_idx: key, 
            allregions: cluster.allregions,
            region: cluster.region,
            country: cluster.country,
            searchtext: cluster.searchtext,
            label1: cluster['lineage'],
            count: coldates.length,
            varcount: cluster['sampled_variants'],
            sampled_varcount: labels.filter(x => x.substring(0,9) !== 'unsampled').length,
            first_date,
            last_date,
            pdist: cluster.pdist,
            rdist: cluster.rdist,
            coldate: last_date,
            x: 0,
            x1: 0,
            x2: 0,
            y: 0,
            max_ndiffs: tip_stats.max_ndiffs,
            mean_ndiffs: tip_stats.mean_ndiffs,
            nsamples: tip_stats.nsamples,
            mutations: tip_stats.mutations,
            residual: tip_stats.mean_ndiffs - exp_diffs,
            mcoldate: first_date.addDays(mean_date),
            infections: tip_stats.infections
          });
        } else {
          const labels = Object.keys(cluster['nodes']);
          const root = tip_labels.filter(v => v === cluster['lineage'])[0];
          if (root !== undefined) {
            tips = map_tips(key, labels, root, tips, tip_labels, cluster);
          } else {
            console.log('Failed to match cluster index', key, 'to a tip in the tree');
          }
        }
      }
    } catch (err) {
      console.error('Error processing cluster key:', key, err.message);
    } finally {
      pipeline.resume();
    }
  });

  await new Promise((resolve, reject) => {
    pipeline.on('end', resolve);
    pipeline.on('error', reject);
  });

  // final flushes
  await flushClusters();
  await flushBeaddata();
  await flushFlat();

  console.log(`Created ${totalClusters} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__CLUSTERS}`);
  console.log(`Created ${totalBeaddata} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__BEADDATA}`);
  console.log(`Created ${totalFlat} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__FLAT_DATA}`);

  const lineage_stats = Object.entries(global.dbstats['lineages']).map(([key, value]) => ({ _id: key, ...value }));
  let res = await db.collection($COLLECTION__DBSTATS).insertMany(lineage_stats, { ordered: false });
  console.log(`Created ${(res.insertedCount ?? Object.keys(res.insertedIds || {}).length)} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__DBSTATS}`);
  delete lineage_stats;

  console.log('Writing accn_to_cid map...');
  {
    const entries = Object.entries(accn_to_cid);
    let chunk = [], written = 0;
    for (const [k, v] of entries) {
      const j = {}; j[k] = v;
      chunk.push(j);
      if (chunk.length >= BATCH_MAP_CHUNKS) { written += await safeInsertMany(db.collection($COLLECTION__ACCN_TO_CID), chunk); chunk = []; }
    }
    if (chunk.length) written += await safeInsertMany(db.collection($COLLECTION__ACCN_TO_CID), chunk);
    console.log(`Created ${written} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__ACCN_TO_CID}`);
  }

  console.log('Writing lineage_to_cid map...');
  {
    const entries = Object.entries(lineage_to_cid);
    let chunk = [], written = 0;
    for (const [k, v] of entries) {
      const j = {}; j[k] = v;
      chunk.push(j);
      if (chunk.length >= BATCH_MAP_CHUNKS) { written += await safeInsertMany(db.collection($COLLECTION__LINEAGE_TO_CID), chunk); chunk = []; }
    }
    if (chunk.length) written += await safeInsertMany(db.collection($COLLECTION__LINEAGE_TO_CID), chunk);
    console.log(`Created ${written} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__LINEAGE_TO_CID}`);
  }
  delete lineage_to_cid;

  res = await db.collection($COLLECTION__REGION_MAP).insertOne(global.region_map);
  console.log(`Created ${res.insertedCount ?? 1} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__REGION_MAP}`);
  delete global.region_map;

  res = await db.collection($COLLECTION__DF_TREE).insertMany(global.df, { ordered: false });
  console.log(`Created ${(res.insertedCount ?? Object.keys(res.insertedIds || {}).length)} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__DF_TREE}`);
  delete global.df;

  res = await db.collection($COLLECTION__XBB_TREE).insertMany(global.df_xbb, { ordered: false });
  console.log(`Created ${(res.insertedCount ?? Object.keys(res.insertedIds || {}).length)} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__XBB_TREE}`);
  delete global.df_xbb;

  const autocomplete_data = Object.keys(accn_to_cid).sort()
    .concat(Object.keys(lineage_to_cid).sort())
    .map(accn => [global.normalize(accn), accn]);

  res = await db.collection($COLLECTION__AUTOCOMPLETE_DATA).insertMany(
    autocomplete_data.map(([norm, accn]) => ({ norm, accn })),
    { ordered: false }
  );
  console.log(`Created ${(res.insertedCount ?? Object.keys(res.insertedIds || {}).length)} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__AUTOCOMPLETE_DATA}`);
  delete autocomplete_data;
  delete accn_to_cid;

  res = await db.collection($COLLECTION__TIPS).insertMany(tips, { ordered: false });
  console.log(`Created ${(res.insertedCount ?? Object.keys(res.insertedIds || {}).length)} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__TIPS}`);
  delete tips;

  recombinant_tips_partial.sort((a, b) => a.first_date - b.first_date);
  for (let i = 0; i < recombinant_tips_partial.length; i++) {
    recombinant_tips_partial[i].y = i;
  }
  if (recombinant_tips_partial.length > 0) {
    res = await db.collection($COLLECTION__RECOMBINANT_TIPS).insertMany(recombinant_tips_partial, { ordered: false });
    console.log(`Created ${(res.insertedCount ?? Object.keys(res.insertedIds || {}).length)} documents in ${$ACTIVE_DATABASE}.${$COLLECTION__RECOMBINANT_TIPS}`);
  }
  delete recombinant_tips_partial;

  await client.close();
  console.log('Sucessfully wrote to database and closed connection!');
}

const question_warning = 'This will delete records within database ' + $ACTIVE_DATABASE + '. Are you sure? (Y/n)';
const userConfirmationCallback = function (answer) {
  if (answer === 'Y') {
    rl.close();
    updateDatabase().catch(async (error) => {
      console.log(`ERROR: ${error}`);
      process.exit(1);
    });
  } else if (answer === 'n') {
    rl.close();
    console.log('EXITING');
    process.exit();
  } else {
    rl.question(question_warning, userConfirmationCallback);
  }
};

if ($NODE_ENV === 'PROD' || $NODE_ENV === 'DEV') {
  rl.close();
  updateDatabase().catch(async (error) => {
    console.log(`ERROR: ${error}`);
    process.exit(1);
  });
} else {
  rl.question(question_warning, userConfirmationCallback);
}
