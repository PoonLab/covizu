const fs = require('fs');

// capture the NODE_ENV env var specified on command-line (if any)
// NODE_ENV='TEST' node server.js
var $HTTP_PORT;
var $HTTPS_PORT;
var $DATA_FOLDER;
var $SSL_CREDENTIALS;
var $NODE_ENV = process.env.NODE_ENV;

// based on the NODE_ENV env var we read a specific .env file

if ($NODE_ENV == 'TEST'){
    require('dotenv').config({ path: '.env.test' });
}
else if ($NODE_ENV == 'PROD'){
    require('dotenv').config({ path: '.env.prod' });
}
else{
    require('dotenv').config({ path: '.env.dev' });
}


$NODE_ENV = process.env.NODE_ENV;//in case we didnt specify the NODE_ENV on the command line as well as inside the .env file this catches it
$DATA_FOLDER = process.env.DATA_FOLDER;
$HTTP_PORT = process.env.HTTP_PORT;
$HTTPS_PORT = process.env.HTTPS_PORT;

if ($NODE_ENV=='PROD') {
    try {
        $SSL_CREDENTIALS = {
            key: fs.readFileSync(process.env.PRVTKEY), 
            cert: fs.readFileSync(process.env.CRT)
        }
    } 
    catch (e) {
        console.error("PROD server requires SSL encryption but .env file is missing PRVTKEY or CRT");
        throw new Error(e);
    }
}

if(!$HTTP_PORT){
    console.warn(".env is missing HTTP_PORT. Defaulting to 8001")
    $HTTP_PORT = 8001;
}

if(!$HTTPS_PORT && $NODE_ENV=='PROD'){
    console.warn(".env is missing HTTPS_PORT. Defaulting to 8002")
    $HTTPS_PORT = 8002;
}
if(!$DATA_FOLDER){
    console.warn('.env is missing DATA_FOLDER env variable. Defaulting to data/')
    $DATA_FOLDER = 'data'
}
if(!$NODE_ENV){
    console.warn('.env is missing NODE_ENV. Defaulting to DEV')
    $NODE_ENV = 'DEV';
}


module.exports = {
    $HTTP_PORT,
    $HTTPS_PORT,
    $DATA_FOLDER,
    $NODE_ENV,
    $SSL_CREDENTIALS
}
