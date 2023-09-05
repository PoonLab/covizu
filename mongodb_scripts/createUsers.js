
/**
A nodeJS script to create the user accounts for the MongoDB database
You must disable authentication for the mongodb server before running this script.
To disable authentication, edit the last two lines of the mongod.conf file
The MongoDB config files are usually located at
MAC INTEL /usr/local/etc/mongod.conf
MAC M1    /opt/homebrew/etc/mongod.conf

Add/Remove the following two lines at the bottom of the config file to enable/disable authentication
    security:
    authorization: enabled
Use a # to comment out the line

After disabling authentication, restart the mongodb service
    brew services restart mongodb/brew/mongodb-community
    sudo systemctl restart mongodb
*/
const MongoClient = require('mongodb').MongoClient;

const {
    $ADMIN_USERNAME,
    $ADMIN_PASSWORD,
    $COVIZU_USERNAME,
    $COVIZU_PASSWORD,
    $DB_URL,
    $DATABASE__PRIMARY,
    $DATABASE__SECONDARY
} = require('../config/dbconfig')

const url = `mongodb://${$DB_URL}`;

// Create a new user account
async function createUserAccounts() {
    const client = new MongoClient(url);

    try {
        // Connect to MongoDB server
        await client.connect();

        // Access the admin database
        const adminDb = client.db('admin');

        let result;

        result = await adminDb.command({
            createUser: $ADMIN_USERNAME,
            pwd: $ADMIN_PASSWORD,
            roles: [
                { role: "root", db: $ADMIN_USERNAME },
                { role: "dbAdmin", db: $DATABASE__PRIMARY },
                { role: "dbAdmin", db: $DATABASE__SECONDARY },
            ]
        });
        console.log('Admin user account created successfully:', result);


        result = await adminDb.command({
            createUser: $COVIZU_USERNAME,
            pwd: $COVIZU_PASSWORD,
            roles: [
                { role: "read", db: $DATABASE__PRIMARY },
                { role: "read", db: $DATABASE__SECONDARY },
            ]
        });
        console.log('User covizu account created successfully:', result);

    } catch (error) {
        console.error('Error creating user account:', error);
    } finally {
        // Close the connection
        client.close();
    }
}

// Call the function to create the user account
createUserAccounts();
