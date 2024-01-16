const { parse_clusters, map_clusters_to_tips, index_accessions, index_lineage } = require('./server/parseCluster')
const { readTree } = require('./server/phylo')
const fs = require('fs');

const { $DATA_FOLDER } = require('./config/config');
const dbstats = require(`./${$DATA_FOLDER}/dbstats.json`)
require("./globalVariables")

class TestDBManager {

    constructor() {/**Nothing to initialize in the TestDBManager */}

    get_edgeList(cindex) {
        return new Promise((resolve,reject)=>{
            resolve(global.beaddata[cindex].edgelist)
        })
    }

    get_points(cindex) {
        return new Promise((resolve,reject)=>{
            resolve(global.beaddata[cindex].points);
        })
    }

    get_variants(cindex) {
        return new Promise((resolve,reject)=>{
            resolve(global.beaddata[cindex].variants);
        })
    }

    get_tips() {
        return new Promise((resolve,reject)=>{
            resolve(global.tips);
        });
    }

    get_df() {
        return new Promise((resolve,reject)=>{
            resolve(global.df);
        });
    }

    get_xbb_df() {
        return new Promise((resolve,reject)=>{
            resolve(global.df_xbb);
        });
    }

    get_regionMap() {
        return new Promise((resolve,reject)=>{
            resolve(global.region_map);
        })
    }

    get_lineage(cindex) {
        return new Promise((resolve,reject)=>{
            resolve(global.clusters[cindex].lineage);
        })
    }

    get_display(lineage) {
        return new Promise((resolve,reject)=>{
            var rawLineage = dbstats["lineages"][lineage]["raw_lineage"];
            if (rawLineage.startsWith("XBB"))
                resolve(["XBB Lineages"]);
            else if (rawLineage.startsWith("X"))
                resolve(["Other Recombinants"]);
            else
                resolve(["Non-Recombinants"]);
        });
    }

    get_accession(accession) {
        return new Promise((resolve,reject)=>{
            resolve(global.accn_to_cid[accession]);
        })
    }

    get_cid() {
        return new Promise((resolve,reject)=>{
            resolve (global.accn_to_cid)
        })
    }

    get_lineageToCid() {
        return new Promise((resolve,reject)=>{
            resolve (global.lineage_to_cid)
        })
    }

    // is global.recombinant_tips an empty array in the test-case?
    get_recombinantTips() {
        return new Promise((resolve,reject)=>{
            resolve(global.recombinant_tips)
        })
    }

    get_searchHits(start_date, end_date, query) {
        return new Promise((resolve,reject)=>{

            // Flatten the json data to an array with bead data only
            let flat_data = global.beaddata.map(bead => bead.points).flat();
            //Find all the beads that are a hit. Convert text_query to lower case and checks to see if there is a match
            let search_hits = flat_data.filter(
                (bead) => 
                (
                    bead.accessions.some(accession => accession.toLowerCase().includes(query)) || 
                    bead.labels.some(label => label.toLowerCase().includes(query))
                    ) && 
                    (bead.x >= start_date) && 
                    (bead.x <= end_date)
                );
            resolve(search_hits);
        })
    }

    get_hits(term) {
        return new Promise((resolve,reject)=>{
            if (!/\d/.test(term)) {
                if (global.prefix.test(term)) {
                    resolve(global.autocomplete_data.slice(0, global.MIN_RESULTS).map(e=>e[1]));
                } else {
                    resolve([]);
                }
            } else {
                const result = global.autocomplete_data.filter(array => array[0].indexOf(normalize(term)) > -1);
                resolve(result.slice(0, global.MIN_RESULTS).map(e=>e[1]));
            }
        })

    }

    load_globalData() {
        global.clusters = require(`./${$DATA_FOLDER}/clusters.json`);
        global.beaddata = parse_clusters(global.clusters)
        try {
            global.tree = fs.readFileSync(`./${$DATA_FOLDER}/timetree.nwk`, 'utf8');
        }
        catch (e) {
            console.log('Error:', e.stack);
        }
        try {
            global.xbbtree = fs.readFileSync(`./${$DATA_FOLDER}/xbbtree.nwk`, 'utf8');
        }
        catch (e) {
            console.log('Error:', e.stack);
        }
        global.df = readTree(global.tree)
        global.df_xbb = readTree(global.xbbtree)
        const { tips, tips_xbb, recombinant_tips } = map_clusters_to_tips(global.df, df_xbb, global.clusters);
        global.tips = tips;
        global.recombinant_tips = recombinant_tips;
        global.accn_to_cid = index_accessions(global.clusters);
        global.lineage_to_cid = index_lineage(global.clusters);
        global.autocomplete_data = Object.keys(global.accn_to_cid).sort()
        .concat(Object.keys(global.lineage_to_cid).sort())
        .map(accn => [global.normalize(accn), accn]);
    }
}

module.exports = { TestDBManager };