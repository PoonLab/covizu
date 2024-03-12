/**
A nodeJS script to delete all exisiting user accounts in the MongoDB database
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
    $COVIZU_USERNAME,
    $DB_URL,
} = require('../config/dbconfig')

const url = `mongodb://${$DB_URL}`;

// Find all user accounts
async function deleteAllUserAccounts() {
    const client = new MongoClient(url);
    let existingUsers = [];
    // Connect to MongoDB server
    await client.connect();
    // Access the admin database
    const adminDb = client.db($ADMIN_USERNAME);

    // Get the list of user accounts
    const users = await adminDb.command({ usersInfo: 1 });

    console.log('User accounts:', users.users);

    try {
        await adminDb.command({ dropUser: $ADMIN_USERNAME});
        console.log(`Deleted user account ${$ADMIN_USERNAME}.${$ADMIN_USERNAME}`)
    }
    catch(error){
        console.log(`Failed to delete user account ${$ADMIN_USERNAME}: ${error}`);
    }
    try
    {
        await adminDb.command({ dropUser: $COVIZU_USERNAME});
        console.log(`Deleted user account ${$ADMIN_USERNAME}.${$COVIZU_USERNAME}`)
    }
    catch(error){
        console.log(`Failed to delete user account ${$COVIZU_USERNAME}: ${error}`);
    }

    await client.close();
}
deleteAllUserAccounts();