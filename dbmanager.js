const { readTree } = require('./server/phylo')
const fs = require('fs');
const { MongoClient } = require("mongodb");

const { $DATA_FOLDER } = require('./config/config');
const {
    $COLLECTION__CLUSTERS,
    $COLLECTION__BEADDATA,
    $COLLECTION__TIPS,
    $COLLECTION__RECOMBINANT_TIPS,
    $COLLECTION__ACCN_TO_CID,
    $COLLECTION__LINEAGE_TO_CID,
    $COLLECTION__REGION_MAP,
    $COLLECTION__DF_TREE,
    $COLLECTION__AUTOCOMPLETE_DATA,
    $COLLECTION__FLAT_DATA
} = require("./config/dbconfig")

require("./globalVariables")

class DBManager {
    url_;
    dbName_;
    db_;
    connection_;

    constructor() {
        this.url_ = null;
        this.dbName_ = null;
        this.connection_ = null;
        this.db_ = null;
    }

    set_dbUrl(url) {
        this.url_ = url;
        this.connection_ = new MongoClient(this.url_);
    }

    set_dbName(dbName) {
        this.dbName_ = dbName;
    }
    get_edgeList(cindex) {
        console.log("GETTING EDGELIST ", cindex)
        cindex = +cindex;
        return this.db_.collection($COLLECTION__BEADDATA).find().skip(cindex).limit(1).toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_edgelist error retrieving edgelist for cindex=${cindex}`,err);
                return;
            }
            if(!docs || !docs.length)
            {
                console.log(`dbmanager::get_edgelist found no records for cindex=${cindex}`);
                return [];
            }
            console.log(docs)
            return docs[0].edgelist;
        })
    }
    get_points(cindex) {
        cindex = +cindex;
        return this.db_.collection($COLLECTION__BEADDATA).find().skip(cindex).limit(1).toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_points error retrieving points for cindex=${cindex}`,err);
                return;
            }
            if(!docs || !docs.length)
            {
                console.log(`dbmanager::get_points found no records for cindex=${cindex}`);
                return [];
            }
            console.log(docs)
            return docs[0].points;
        })
    }

    get_variants(cindex) {
        cindex = +cindex;
        return this.db_.collection($COLLECTION__BEADDATA).find().skip(cindex).limit(1).toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::variants error retrieving points for cindex=${cindex}`,err);
                return;
            }
            if(!docs || !docs.length)
            {
                console.log(`dbmanager::variants found no records for cindex=${cindex}`);
                return [];
            }
            console.log(docs)
            return docs[0].variants;
        })
    }

    get_tips() {
        console.log("GETTING TIPS")
        return this.db_.collection($COLLECTION__TIPS).find().toArray()
    }

    get_df() {
        return this.db_.collection($COLLECTION__DF_TREE).find().toArray();
    }

    get_regionMap() {
        return this.db_.collection($COLLECTION__REGION_MAP).find().toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_regionMap error`,err);
                return ;
            }
            return docs[0];
        });
        // return global.region_map;
    }

    get_lineage(cindex) {
        console.log("GETTING EDGELIST ", cindex)
        cindex = +cindex;
        return this.db_.collection($COLLECTION__CLUSTERS).find().skip(cindex).limit(1).toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_lineage error retrieving edgelist for cindex=${cindex}`,err);
                return;
            }
            if(!docs || !docs.length)
            {
                console.log(`dbmanager::get_lineage found no records for cindex=${cindex}`);
                return [];
            }
            console.log(docs)
            return docs[0].lineage;
        })
    }

    get_accession(accession) {
        const query = { [accession]: { $exists: true } };
        return this.db_.collection($COLLECTION__ACCN_TO_CID)
        .findOne(query).then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_accession - error accession=${accession}`);
                return;
            }
            if(docs)
                return docs[accession];
            else{
                console.log(`dbmanager::get_accession - no documents found matching accession=${accession}`);
                return;
            }
        })
        // return global.accn_to_cid[accession]
    }

    get_cid() {
        return this.db_.collection($COLLECTION__ACCN_TO_CID)
        .aggregate([ 
            // {$limit: 50}, // in theory, this endpoint works. but, do we want to serve a 60MB file?
            { $project: { _id: 0 } },
            { $group: { _id: null, mergedObject: { $mergeObjects: "$$ROOT" } } }, 
            { $replaceRoot: { newRoot: "$mergedObject" } },
        ])
        .toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_cid error`,err);
                return ;
            }
            if(docs && docs.length){
                console.log("found docs",docs);
                return docs;
            }
            else{
                console.log(`dbmanager::get_cid no documents found`)
                return;
            }
        })
        // return (global.accn_to_cid)
    }

    get_lineageToCid() {
        return this.db_.collection($COLLECTION__LINEAGE_TO_CID)
        .aggregate([ 
            { $project: { _id: 0 } },
            { $group: { _id: null, mergedObject: { $mergeObjects: "$$ROOT" } } }, 
            { $replaceRoot: { newRoot: "$mergedObject" } },
        ])
        .toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_lineageToCid error`,err);
                return ;
            }
            if(docs && docs.length == 1){
                // console.log("found docs",docs);
                return docs[0];
            }
            else{
                console.log(`dbmanager::get_lineageToCid no documents found`);
                return ;
            }
        })
        // return global.lineage_to_cid
    }

    get_recombinantTips() {
        return this.db_.collection($COLLECTION__RECOMBINANT_TIPS).find().toArray().then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_recombinantTips error`,err)
                return;
            }
            return docs;

        })
        // return global.recombinant_tips
    }

    // get_searchHits(start_date, end_date, query) {
    //     // Flatten the json data to an array with bead data only
    //     let flat_data = global.beaddata.map(bead => bead.points).flat();

    //     //Find all the beads that are a hit. Convert text_query to lower case and checks to see if there is a match
    //     let search_hits = flat_data.filter(function (bead) {
    //         let temp = (bead.accessions.some(accession => (accession.toLowerCase()).includes(query)) ||
    //             bead.labels.some(label => (label.toLowerCase()).includes(query))) &&
    //             (bead.x >= start_date && bead.x <= end_date);
    //         return temp;
    //     });
    //     return search_hits;
    // }

    get_searchHits(start_date,end_date,term){
        const regexTerm = new RegExp(term,'i');
        const query = {
            $or: [
                { labels: { $regex: regexTerm } },
                { accessions: { $regex: regexTerm } }
            ],
            start_date: { $gt: new Date(start_date) },
            end_date: { $lt: new Date(end_date) }
        }
        return this.db_.collection($COLLECTION__FLAT_DATA).find(query).toArray()
        .then((docs,err)=>{
            if(err){
                console.log(`dbmanager::get_searchHits error ${start_date},${end_date},${query}`,err);
                return [];
            }
            if(!docs || !docs.length){
                console.log(`dbmanager::get_searchHits no docs found ${start_date},${end_date},${query}`);
                return [];
            }
            // console.log("FOUND get_searchHIts", docs);
            return docs;
        })
    }
    
    get_hits_callback(docs,err,term){
        if(err){
            console.log(`dbmanager::get_hits error ${term}`,err);
            return [];
        }
        if(!docs || !docs.length){
            console.log(`dbmanager::get_hits found 0 docs ${term}`);
            return [];
        }
        // console.log(docs.map(e=>e['accn']));
        return docs.map(e=>e['accn']);
    }

    get_hits(term) {
        if (!/\d/.test(term)) {
            if (global.prefix.test(term)) {
                return this.db_.collection($COLLECTION__AUTOCOMPLETE_DATA).find().limit(global.MIN_RESULTS).toArray()
                .then((docs,err)=>{return this.get_hits_callback(docs,err,term)})
            } else {
                return([]);
            }
        } else {
            const query = { "norm": { $regex: term, $options: "i" } };
            return this.db_.collection($COLLECTION__AUTOCOMPLETE_DATA).find(query).limit(global.MIN_RESULTS).toArray()
            .then((docs,err)=>{return this.get_hits_callback(docs,err,term)})
        }
    }

    async load_globalData() {
        // console.log("Loading clusters")
        // this.db_.collection($COLLECTION__CLUSTERS).find().toArray().then(docs=>{
        //     global.clusters = docs;
        //     console.log("Loading beaddata")
        //     this.db_.collection($COLLECTION__BEADDATA).find().toArray().then(docs=>{
        //         global.beaddata = docs;
        //         console.log("Loading tips")
        //         this.db_.collection($COLLECTION__TIPS).find().toArray().then(docs=>{
        //             global.tips =  docs;
        //             console.log("Loading recombinant_tips");
        //             this.db_.collection($COLLECTION__RECOMBINANT_TIPS).find().toArray().then(docs=>{
        //                 global.recombinant_tips = docs;
        //                 console.log("Loading accn_to_cid")
        //                 this.db_.collection($COLLECTION__ACCN_TO_CID).find().toArray().then(docs=>{
        //                     global.accn_to_cid = docs;
        //                     console.log("Loading lineage_to_cid")
        //                     this.db_.collection($COLLECTION__LINEAGE_TO_CID).find().toArray().then(docs=>{
        //                         global.lineage_to_cid = docs;
        //                         console.log("Preparing global.autocomplete_data")
        //                         global.autocomplete_data = Object.keys(global.accn_to_cid).sort().concat(Object.keys(global.lineage_to_cid).sort()).map(accn => [
        //                             global.normalize(accn), accn
        //                         ]);
        //                     })
        //                 })
        //             })
        //         });
        //     })
        // });
        

        // console.log("Loading clusters")
        // global.clusters = await this.db_.collection($COLLECTION__CLUSTERS).find().toArray();

        // console.log("Loading beaddata")
        // global.beaddata = await this.db_.collection($COLLECTION__BEADDATA).find().toArray();
        
        // console.log("Loading tips")
        // global.tips = await this.db_.collection($COLLECTION__TIPS).find().toArray();

        // console.log("Loading recombinant_tips")
        // global.recombinant_tips = await this.db_.collection($COLLECTION__RECOMBINANT_TIPS).find().toArray();

        // console.log("Loading accn_to_cid")
        // global.accn_to_cid = await this.db_.collection($COLLECTION__ACCN_TO_CID).find().toArray();

        // console.log("Loading lineage_to_cid")
        // global.lineage_to_cid = await this.db_.collection($COLLECTION__LINEAGE_TO_CID).find().toArray();

        // await this.collClusters_.find().toArray()
        //     .then((docs) => {
        //         global.clusters = docs;
        //         global.beaddata = parse_clusters(docs);
        //     })

        // const { tips, recombinant_tips } = map_clusters_to_tips(global.df, global.clusters);
        // global.tips = tips;
        // global.recombinant_tips = recombinant_tips;
        // global.accn_to_cid = index_accessions(global.clusters);
        // global.lineage_to_cid = index_lineage(global.clusters);
        // global.autocomplete_data = Object.keys(global.accn_to_cid).sort().concat(Object.keys(global.lineage_to_cid).sort()).map(accn => [
        //     global.normalize(accn), accn
        // ]);
    }

    async start() {
        try {
            if (this.connection_)
                await this.connection_.close();
            if (!this.url_)
                throw new Error("Missing dbUrl")
            if (!this.dbName_)
                throw new Error("Missing dbName")
            await this.connection_.connect();
            this.db_ = this.connection_.db(this.dbName_);
            console.log("Connected to database");
        } catch (error) {
            console.log("DbManager Error: ", error);
        }
    }

    async close() {
        if (!this.connection_) return;
        try {
            await this.connection_.close();
            console.log('Disconnected from the database');
        } catch (error) {
            console.error('Error disconnecting from the database:', error);
        }
    }
}

module.exports = { DBManager };
