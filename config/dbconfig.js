const path = require('path');

const $PROJECT_ROOT = path.resolve(__dirname+"/../");
const $JSON_DATA_FOLDER = "data"
const $DATABASE__PRIMARY = "covizu_1";
const $DATABASE__SECONDARY = "covizu_2";
const $COLLECTION__CLUSTERS = "clusters";
const $COLLECTION__BEADDATA = "beaddata";
const $COLLECTION__TIPS = "tips";
const $COLLECTION__RECOMBINANT_TIPS = "recombinant_tips";
const $COLLECTION__ACCN_TO_CID = "accn_to_cid";
const $COLLECTION__LINEAGE_TO_CID = "lineage_to_cid";
const $COLLECTION__REGION_MAP = "region_map";
const $COLLECTION__DF_TREE = "df_tree";
const $COLLECTION__AUTOCOMPLETE_DATA = "autocomplete_data";
const $COLLECTION__FLAT_DATA = "flat_data";

const $JSONFILE__CLUSTERS = "clusters.json";
const $NWKFILE__TREE = "timetree.nwk";

var $NODE_ENV = process.env.NODE_ENV;

if ($NODE_ENV == 'TEST'){
    require('dotenv').config({ path: '.env.test' });
}
else if ($NODE_ENV == 'PROD'){
    require('dotenv').config({ path: '.env.prod' });
}
else{
    require('dotenv').config({ path: '.env.dev' });
}

let $ACTIVE_DATABASE;
if($NODE_ENV != 'TEST'){
    const $DBNUMBER = process.env.DBNUMBER;
    if($DBNUMBER == 1)
    {
        console.log("Setting active database to ", $DATABASE__PRIMARY);
        $ACTIVE_DATABASE = $DATABASE__PRIMARY;
    }
    else if($DBNUMBER == 2)
    {
        console.log("Setting active database to ", $DATABASE__SECONDARY);
        $ACTIVE_DATABASE = $DATABASE__SECONDARY;
    }
    else
    {
        $ACTIVE_DATABASE = undefined;
    }
}

const $ADMIN_USERNAME = encodeURIComponent(process.env.ADMIN_USERNAME);
const $ADMIN_PASSWORD = encodeURIComponent(process.env.ADMIN_PASSWORD);
const $COVIZU_USERNAME = encodeURIComponent(process.env.COVIZU_USERNAME);
const $COVIZU_PASSWORD = encodeURIComponent(process.env.COVIZU_PASSWORD);
const $AUTHMECHANISM = "DEFAULT";
const $DB_URL = process.env.DB_URL;

// this covizu_admin connection is required to create the database/collections/documents.
const $ADMIN_CONNECTION_URI = `mongodb://${$ADMIN_USERNAME}:${$ADMIN_PASSWORD}@${$DB_URL}/?authMechanism=${$AUTHMECHANISM}`;
// this covizu connection will be used for reading data from the database and delivering to the webserver and the frontend
const $COVIZU_CONNECTION_URI = `mongodb://${$COVIZU_USERNAME}:${$COVIZU_PASSWORD}@${$DB_URL}/?authMechanism=${$AUTHMECHANISM}`;

module.exports = {
    $DATABASE__PRIMARY,
    $DATABASE__SECONDARY,
    $ADMIN_USERNAME,
    $ADMIN_PASSWORD,
    $COVIZU_USERNAME,
    $COVIZU_PASSWORD,
    $DB_URL,
    $ADMIN_CONNECTION_URI,
    $COVIZU_CONNECTION_URI,
    $ACTIVE_DATABASE,
    $COLLECTION__CLUSTERS,
    $COLLECTION__BEADDATA,
    $COLLECTION__TIPS,
    $COLLECTION__RECOMBINANT_TIPS,
    $COLLECTION__ACCN_TO_CID,
    $COLLECTION__LINEAGE_TO_CID,
    $COLLECTION__REGION_MAP,
    $COLLECTION__DF_TREE,
    $COLLECTION__AUTOCOMPLETE_DATA,
    $COLLECTION__FLAT_DATA,
    $PROJECT_ROOT,
    $JSON_DATA_FOLDER,
    $JSONFILE__CLUSTERS,
    $NWKFILE__TREE,
}