// capture the NODE_ENV env var specified on command-line (if any)
// NODE_ENV='TEST' node server.js
var $HTTP_PORT;
var $DATA_FOLDER;
var $NODE_ENV = process.env.NODE_ENV;

// based on the NODE_ENV env var we read a specific .env file
if ($NODE_ENV == 'TEST') {
    require('dotenv').config({ path: '.env.test' });
}
else if ($NODE_ENV == 'PROD') {
    require('dotenv').config({ path: '.env.prod' });
}
else {
    require('dotenv').config({ path: '.env.dev' });
}

$NODE_ENV = process.env.NODE_ENV;//in case we didnt specify the NODE_ENV on the command line...
$DATA_FOLDER = process.env.DATA_FOLDER;
$HTTP_PORT = process.env.HTTP_PORT;

if (!$HTTP_PORT) {
    console.warn(".env is missing HTTP_PORT. Defaulting to 8001")
    $HTTP_PORT = 8001;
}

if (!$DATA_FOLDER) {
    console.warn('.env is missing DATA_FOLDER env variable. Defaulting to data/')
    $DATA_FOLDER = 'data'
}
if (!$NODE_ENV) {
    console.warn('.env is missing NODE_ENV. Defaulting to DEV')
    $NODE_ENV = 'DEV';
}

module.exports = {
    $HTTP_PORT,
    $DATA_FOLDER,
    $NODE_ENV,
}
