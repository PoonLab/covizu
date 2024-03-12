## Installing `MongoDB Community Version`
### On Ubuntu
  * [Full Instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/)
### On MacOS
  * [Full Instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-os-x/)
### On Windows
  * [Full Instructions](https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-windows/)
________________

## Directories and Files (`config`, `logs`, `data`)
These are the _default locations_ for some important `mongodb` files. 
### On Ubuntu
  * __config file__ : `/etc/mongod.conf`
  * __log directory__ : `/var/log/mongodb`
  * __data directotry__ : `/var/lib/mongodb`
### On MacOS Intel Processor
  * __config file__ : `/usr/local/etc/mongod.conf`
  * __log directory__ : `/usr/local/var/log/mongodb`
  * __data directotry__ : `/usr/local/var/mongodb`
### On MacOS M1 Processor
  * __config file__ :   `/opt/homebrew/etc/mongod.conf`
  * __log directory__ :   `/opt/homebrew/var/log/mongodb`
  * __data directory__ :   `/opt/homebrew/var/mongodb`
________________

## To start/stop/restart the `MongoDB  daemon`
_When running as a service, the default config file will be used_ 
### On Ubuntu
  * `$ sudo systemctl start mongod`
  * `$ sudo systemctl stop mongod`
  * `$ sudo systemctl restart mongod`

### On MacOS
  * `$ brew services start mongodb/brew/mongodb-community`
  * `$ brew services stop mongodb/brew/mongodb-community`
  * `$ brew services restart mongodb/brew/mongodb-community`
________________

## To manually start `MongoDB` with a custom config file
### For MacOS Intel Processors
* `$ mongod --config /usr/local/etc/mongod.conf --fork`
### For MacOS M1 processors
* `$ mongod --config /opt/homebrew/etc/mongod.conf --fork`
* _NOTE_: This could be useful when needing to run a `test` database with separate `log` files and `data` directory. Or we could use the same `log` files and the same `data` directory and instead just create a new database `covizu_test` within the Mongo application.
________________

## To open the Command Line Tool for `MongoDB` a.k.a `mongosh`
* `$ mongosh`

## Database Authorization
* By default MongoDB does not seem to enable any authorization. Meaning anyone can connect to our database on `port:27017` and CRUD our data.
* However, it might be different for you. When trying to open the `cli` if you get an error that looks like this
```
Current Mongosh Log ID:	63fc56344c05787dbf853d87
Connecting to:mongodb://127.0.0.1:27017/?directConnection=true&serverSelectionTimeoutMS=2000&appName=mongosh+1.7.1
MongoNetworkError: connect ECONNREFUSED 127.0.0.1:27017
```
then it probably means your `MongoDB` has locked you out. You need to disable authorization. 
* To enable/disable authorization, open the afore-mentioned `mongod.conf` configuraiton file.
* Add (remove) the following lines to the bottom of the config file to enable (disable) authorization
```
security:
  authorization: enabled
```
#### _NOTE_ : That `whitespace` between `authorization:` and `enabled` might be crucial. 
* After editing those lines, restart the `MongoDB daemon`. 
* To confirm that you've succesfully enabled(disabled) authorization, open `mongosh` without any auth token.
    * `$ mongosh`
* Switch to the `admin` database 
    * `$ use admin`
* View all the databases 
    * `$ show databases`
* If `MongoDB` prints out something like this then you've succesfully disabled (unsuccesfully enabled) authentication.
```
admin> show databases
admin     180.00 KiB
config    108.00 KiB
local     120.00 KiB
```
* If instead, `MongoDB` prints out something like this then you've succesfully enabled (unsuccesfully disabled) authentication.
```
admin> show databases
MongoServerError: command listDatabases requires authentication
```
___________________


## Creating user accounts
* On first install, we need to create a `admin` account, and a `covizu` acccount (readonly). To do so, first disable authorization in the `mongod.conf` configuration file
* Then ensure that your `.env.dev` or `.env.prod` file contains the environment variables for the following
```
DATA_FOLDER='data'
ADMIN_USERNAME='admin'
ADMIN_PASSWORD='supersecretpassword'
COVIZU_USERNAME='covizu'
COVIZU_PASSWORD='supersecretpassword'
DB_URL='localhost:27017'
```
* Before creating the user accounts, first we check for any existing user accounts of the same name as our `ADMIN_USERNAME` or `COVIZU_USERNAME`. Run the following `nodejs` script.
    * `$ npm run delete-users`
* To confirm that there are no user accounts with conflicting names open the `mongosh` 
    * `$ mongosh`
    * `test> use admin`
    * `test> show users`
* After deleting any existing user accounts, run the following script
  * `$ npm run create-users`


* Once your user accounts are created, edit the `config` file to enable authorization and restart the service. 

## Populating the database
* If you haven't already, obtain the following data files from our main server at https://filogeneti.ca/covizu/data/:
* `timetree.nwk`
* `dbstats.json`
* `clusters.json`
and save your local copies under `covizu/data/`.
* To import these data files into our `MongoDB` database run the following script
    * `$ npm run update-db1`
* This will create our primary database and import the `JSON`/`text` records into it.
* Once the import script has finished executing, start the server by executing the following script
    * `$ npm run start-db1`
* Navgiate your browser to `localhost:8001` to verify that it works. 
* If that went well, then run the next script to setup our secondary database
    * `$ npm run update-db2`
* After executing the import script, start the server with a connection to the secondary database
    * `$ npm run start-db2`

## Updating the database
* To update the database, first obtain the new data files from `https://filogeneti.ca/covizu/data` and replace the files in `covizu/data`
* If your `node` server is currently using the primary (secondary) database then update the secondary (primary) database
    * `$ npm run update-db2` (`$ npm run update-db1`) 
* After updating the secondary (primary) database, close the currently running server which is still serving data from the primary(secondary) database and restart the server with a connection to the secondary (primary) database
    * `$ npm run start-db2` (`$ npm run start-db1`)
* On the next batch of new data, update the primary (secondaary) database and restart the node server with a connection to the primary (secondary) database


## Some useful `mongosh` shell commands
* To open `mongo` CLI in non-authenticated mode
    * `$ mongosh`
* To open `mongo` CLI as an `admin`
    * `$ mongosh --username admin --authenticationDatabase admin`

* To view all users
    * Connect as `admin` ()
    * `use admin`
    * `show users`

* To delete an existing user `"myusername"`
    * Connect as `admin`
    * `use admin`
    * `show users` 
    * `db.dropUser('myusername')`

* To view all collections within a database `"mydatabase"` 
    * Connect as `admin`
    * `use mydatabase`
    * `show collections`

* To drop the database `"mydatabase"`
    * Connect as `admin`
    * `use mydatabase`
    * `db.dropDatabase()`

* To drop the collection `"mycollection"` with the database `"mydatabase"`
    * Connect as `admin`
    * `use mydatabase`
    * `db.mycollection.drop()`
