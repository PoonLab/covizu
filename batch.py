"""
what batch.py does
"""
import argparse
import os
import sys
import json
from datetime import datetime
from tempfile import NamedTemporaryFile
import subprocess
from Bio import Phylo

import covizu
from covizu.utils import gisaid_utils
from covizu.utils.progress_utils import Callback
from covizu.utils.batch_utils import (unpack_records, build_timetree, parse_alias,
                                      make_beadplots, get_mutations)
from covizu.utils.seq_utils import SC2Locator


def parse_args():
    """
    parse inputted arguments
    """
    parser = argparse.ArgumentParser(
        description="CoVizu analysis pipeline automation")

    parser.add_argument(
        '--url',
        type=str,
        default=os.environ.get(
            "GISAID_URL",
            None),
        help="URL to download provision file, defaults to environment variable.")
    parser.add_argument(
        '--user',
        type=str,
        default=os.environ.get(
            "GISAID_USER",
            None),
        help="GISAID username, defaults to environment variable.")
    parser.add_argument(
        '--password',
        type=str,
        default=os.environ.get(
            "GISAID_PSWD",
            None),
        help="GISAID password, defaults to environment variable.")

    parser.add_argument(
        "--infile",
        type=str,
        default=None,
        help="input, path to xz-compressed JSON; if not specified, "
        "download xz file from GISAID provision feed.")
    parser.add_argument("--outdir", type=str, default='data/',
                        help="option, path to write output files")

    parser.add_argument(
        '--minlen',
        type=int,
        default=29000,
        help='option, minimum genome length (nt)')
    parser.add_argument(
        '--mindate',
        type=str,
        default='2019-12-01',
        help='option, earliest possible sample collection date (ISO format, default '
        '2019-12-01')
    parser.add_argument(
        '--poisson-cutoff',
        type=float,
        default=0.001,
        help='option, filtering outlying genomes whose distance exceeds the upper '
        'quantile of Poisson distribution (molecular clock).  Default 0.001 '
        'corresponds to 99.9%% cutoff.')

    parser.add_argument(
        '--batchsize',
        type=int,
        default=2000,
        help='option, number of records to batch process with minimap2')
    parser.add_argument(
        '--max-variants',
        type=int,
        default=5000,
        help='option, limit number of variants per lineage (default 5000)')

    parser.add_argument(
        "--ref",
        type=str,
        default=os.path.join(
            covizu.__path__[0],
            "data/NC_045512.fa"),
        help="option, path to FASTA file with reference genome")
    parser.add_argument('--mmbin', type=str, default='minimap2',
                        help="option, path to minimap2 binary executable")
    parser.add_argument('-mmt', "--mmthreads", type=int, default=16,
                        help="option, number of threads for minimap2.")

    parser.add_argument(
        '--misstol',
        type=int,
        default=300,
        help="option, maximum tolerated number of missing bases per "
        "genome (default 300).")
    parser.add_argument(
        "--vcf",
        type=str,
        default=os.path.join(
            covizu.__path__[0],
            "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf"),
        help="Path to VCF file of problematic sites in SARS-COV-2 genome. "
        "Source: https://github.com/W-L/ProblematicSites_SARS-CoV2")

    parser.add_argument('--ft2bin', default='fasttree2',
                        help='option, path to fasttree2 binary executable')

    parser.add_argument(
        '--alias',
        type=str,
        default=os.path.join(
            covizu.__path__[0],
            "data/pango-designation/pango_designation/alias_key.json"),
        help="optional, path to JSON file containing alias.")
    parser.add_argument(
        '--lineages',
        type=str,
        default=os.path.join(
            covizu.__path__[0],
            "data/pango-designation/lineages.csv"),
        help="optional, path to CSV file containing Pango lineage designations.")
    parser.add_argument('--ttbin', default='treetime',
                        help='option, path to treetime binary executable')
    parser.add_argument('--clock', type=float,
                        help='option, specify molecular clock rate for '
                             'constraining Treetime analysis.')
    parser.add_argument(
        '--earliest',
        action='store_true',
        help='option, use earliest sample per lineage for time-scaled '
        'tree; otherwise defaults to most recent samples.')

    parser.add_argument('--datetol', type=float, default=0.1,
                        help='option, exclude tips from time-scaled tree '
                             'with high discordance between estimated and '
                             'known sample collection dates (year units,'
                             'default: 0.1)')

    parser.add_argument('--binpath', type=str, default='rapidnj',
                        help='option, path to RapidNJ binary executable')
    parser.add_argument('--mincount', type=int, default=5000,
                        help='option, minimum number of variants in lineage '
                             'above which MPI processing will be used.')
    parser.add_argument('--machine_file', type=str, default='mfile',
                        help='option, path to machine file for MPI.')
    parser.add_argument("-n", "--nboot", type=int, default=100,
                        help="Number of bootstrap samples, default 100.")

    parser.add_argument(
        "--boot-cutoff",
        type=float,
        default=0.5,
        help="Bootstrap cutoff for consensus tree (default 0.5). "
        "Only used if --cons is specified.")

    parser.add_argument(
        '--use-db',
        action="store_true",
        help="Use a database to store and retrieve features for sequences")
    parser.add_argument(
        '--dbname',
        type=str,
        default=os.environ.get(
            "POSTGRES_DB",
            "gisaid_db"),
        help="Postgresql database name")
    parser.add_argument(
        '--dbhost',
        type=str,
        default=os.environ.get(
            "POSTGRES_HOST",
            "localhost"),
        help="Postgresql database host address")
    parser.add_argument(
        '--dbport',
        type=str,
        default=os.environ.get(
            "POSTGRES_PORT",
            "5432"),
        help="Connection to port number")
    parser.add_argument(
        '--dbuser',
        type=str,
        default=os.environ.get(
            "POSTGRES_USER",
            None),
        help="Postgresl user")
    parser.add_argument(
        '--dbpswd',
        type=str,
        default=os.environ.get(
            "POSTGRES_PASSWORD",
            None),
        help="Postgresl password")

    parser.add_argument("--dry-run", action="store_true",
                        help="Do not upload output files to webserver.")

    return parser.parse_args()


def open_connection(database_connection_parameters):
    """ open connection to database, initialize tables if they don't exist
        :out:
            :cursor: interactive sql object containing tables
    """
    database_conn = psycopg2.connect(**database_connection_parameters)
    connection_cur = database_conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    # create tables if they don't exist
    seqs_table = '''CREATE TABLE IF NOT EXISTS SEQUENCES (accession VARCHAR(255)
                    PRIMARY KEY, qname VARCHAR(255), lineage VARCHAR(255),
                    date VARCHAR(255), location VARCHAR(255),
                    diffs VARCHAR, missing VARCHAR)'''
    connection_cur.execute(seqs_table)

    connection_cur.execute(
        '''CREATE INDEX IF NOT EXISTS accession_index ON SEQUENCES (accession)''')

    cluster_table = '''
    CREATE TABLE IF NOT EXISTS CLUSTERS (
        lineage VARCHAR(255) PRIMARY KEY,
        cluster_data JSONB
    );
    '''
    connection_cur.execute(cluster_table)
    connection_cur.execute(
        '''CREATE INDEX IF NOT EXISTS lineage_index on CLUSTERS (lineage)''')

    updated_lineages_query = '''
    CREATE TABLE IF NOT EXISTS NEW_RECORDS (
        accession VARCHAR(255),
        lineage VARCHAR(255)
    );
    '''
    connection_cur.execute(updated_lineages_query)

    clear_updated_lineages = '''
    DELETE FROM NEW_RECORDS
    '''
    connection_cur.execute(clear_updated_lineages)

    dbstats_table = '''
    CREATE TABLE IF NOT EXISTS DBSTATS (
        lineage VARCHAR(255) PRIMARY KEY,
        dbstats_data JSONB
    );
    '''
    connection_cur.execute(dbstats_table)

    database_conn.commit()
    return connection_cur, database_conn


def process_feed(input_args, input_cur, callback=None):
    """ Process feed data """
    if callback:
        callback("Processing GISAID feed data")
    loader = gisaid_utils.load_gisaid(
        input_args.infile,
        minlen=input_args.minlen,
        mindate=input_args.mindate)
    batcher = gisaid_utils.batch_fasta(loader, input_cur, size=input_args.batchsize)
    aligned = gisaid_utils.extract_features(
        batcher,
        ref_file=input_args.ref,
        cur=input_cur,
        binpath=input_args.mmbin,
        nthread=input_args.mmthreads,
        minlen=input_args.minlen)
    filtered = gisaid_utils.filter_problematic(
        aligned,
        vcf_file=input_args.vcf,
        cutoff=input_args.poisson_cutoff,
        callback=callback)
    return gisaid_utils.sort_by_lineage(filtered, callback=callback)


if __name__ == "__main__":
    args = parse_args()
    cb = Callback()

    CUR = None
    if args.use_db:
        import psycopg2
        import psycopg2.extras
        from psycopg2 import sql
        from psycopg2.errors import DuplicateDatabase

        # Check if database exists
        connection_parameters = {
            "host": args.dbhost,
            "port": args.dbport,
            "user": args.dbuser,
            "password": args.dbpswd,
        }

        CONNECTION = None
        try:
            CONNECTION = psycopg2.connect(**connection_parameters)
            CONNECTION.autocommit = True

            cursor = CONNECTION.cursor()
            cursor.execute(
                sql.SQL('CREATE DATABASE {}').format(
                    sql.Identifier(
                        args.dbname)))
            cb.callback(f"Database {args.dbname} created successfully.")
        except DuplicateDatabase:
            cb.callback(f"Database {args.dbname} already exists.")
        except psycopg2.Error as e:
            cb.callback(f"Error initiating connection to database: {e}")
            sys.exit()
        finally:
            if CONNECTION is not None:
                cursor.close()
                CONNECTION.close()

        connection_parameters['dbname'] = args.dbname
        CUR, conn = open_connection(connection_parameters)

    # check that user has loaded openmpi module
    try:
        subprocess.check_call(['mpirun', '-np', '2', 'ls'],
                              stdout=subprocess.DEVNULL)
    except FileNotFoundError:
        cb.callback(
            "mpirun not loaded - run `module load openmpi/gnu`",
            level='ERROR')
        sys.exit()

    # check that the user has included submodules
    if (
        not os.path.exists(
            os.path.join(
            covizu.__path__[0],
            "data/pango-designation/lineages.csv")) or not os.path.exists(
                os.path.join(
                    covizu.__path__[0],
                    "data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf"))):
        try:
            subprocess.check_call(
                "git submodule init; git submodule update", shell=True)
        except subprocess.CalledProcessError as e:
            cb.callback("Error adding the required submodules")
            sys.exit()

    # update submodules
    try:
        subprocess.check_call(
            "git submodule foreach git pull origin master",
            shell=True)
    except subprocess.CalledProcessError as e:
        cb.callback(f"Could not update submodules, exit {e}", level='ERROR')

    # download xz file if not specified by user
    if args.infile is None:
        cb.callback("No input specified, downloading data from GISAID feed...")
        args.infile = gisaid_utils.download_feed(
            args.url, args.user, args.password)

    # filter data, align genomes, extract features, sort by lineage
    by_lineage = process_feed(args, CUR, cb.callback)

    # separate XBB and other recombinant lineages
    aliases = parse_alias(args.alias)
    designation = {}
    for prefix, truename in aliases.items():
        if isinstance(truename, list):
            designation.update({prefix: {
                'type': 'XBB' if prefix == 'XBB' else 'recombinant',
                'fullname': '/'.join(truename)
            }})
        else:
            designation.update({prefix: {'type': 'XBB' if truename.startswith(
                "XBB") else 'non-recombinant', 'fullname': truename}})

    # use results to partition by_lineage database
    non_recomb = {}
    XBB = {}
    other_recomb = {}
    for lineage, ldata in by_lineage.items():
        # Put unassigned lineages in non-recombinant category
        if lineage.lower() == "unassigned":
            non_recomb.update({lineage: ldata})
            continue

        prefix = lineage.split('.')[0]
        category = designation[prefix]['type']
        if category == 'non-recombinant':
            non_recomb.update({lineage: ldata})
        elif category == 'XBB':
            XBB.update({lineage: ldata})
        else:
            other_recomb.update({lineage: ldata})

    if len(XBB) < 2:
        other_recomb.update(XBB)
        XBB = None  # no point in building a tree

    # reconstruct time-scaled trees
    timetree, residuals = build_timetree(
        non_recomb, args, callback=cb.callback)
    timestamp = datetime.now().isoformat().split('.')[0]
    nwk_file = os.path.join(args.outdir, f'timetree.{timestamp}.nwk')
    with open(nwk_file, 'w', encoding='utf-8') as handle:
        Phylo.write(timetree, file=handle, format='newick')

    xbb_file = os.path.join(args.outdir, f'xbbtree.{timestamp}.nwk')
    # xbb_outgrp.fa is GISAID-derived data and should NOT be committed to repo
    outgrp_file = os.path.join(covizu.__path__[0], "data/xbb_outgrp.fa")
    with open(xbb_file, 'w', encoding='utf-8') as handle:
        if XBB is not None:
            timetree_xbb, residuals_xbb = build_timetree(
                XBB, args, outgroup=outgrp_file, callback=cb.callback
            )
            residuals.update(residuals_xbb)
            Phylo.write(timetree_xbb, file=handle, format='newick')
        # else empty file

    UPDATED_LINEAGES = None
    if args.use_db:
        UPDATED_LINEAGES_QUERY = '''
        SELECT DISTINCT LINEAGE FROM NEW_RECORDS;
        '''
        CUR.execute(UPDATED_LINEAGES_QUERY)
        UPDATED_LINEAGES = [row['lineage'] for row in CUR.fetchall()]

    # clustering analysis of lineages
    result, infection_prediction = make_beadplots(
        by_lineage, args, cb.callback, initial_time=cb.intial_time.timestamp(), updated_lineages=UPDATED_LINEAGES)
    clust_file = os.path.join(args.outdir, f'clusters.{timestamp}.json')
    if args.use_db:
        # Insert all updated records into the database
        for record in result:
            CUR.execute('''
            INSERT INTO CLUSTERS
            VALUES (%s, %s)
            ON CONFLICT (lineage) DO UPDATE
            SET cluster_data = %s
            ''', [record['lineage'], json.dumps(record), json.dumps(record)])

        # Retrieve cluster data for other lineages from the database
        for lineage, _ in by_lineage.items():
            if lineage not in UPDATED_LINEAGES:
                CUR.execute(
                    f"SELECT cluster_data FROM CLUSTERS WHERE lineage = '{lineage}'")
                cluster_info = CUR.fetchone()
                if cluster_info is None:
                    cb.callback(
                        f"Missing CLUSTERS record for lineage {lineage}",
                        level='ERROR')
                    sys.exit()
                result.append(cluster_info['cluster_data'])

    with open(clust_file, 'w', encoding='utf-8') as handle:
        json.dump(result, fp=handle)

    # get mutation info
    locator = SC2Locator()
    mutations = {}
    for lineage, features in get_mutations(by_lineage).items():
        annots = {
            locator.parse_mutation(f): freq for f,
            freq in features.items()}
        mutations.update(
            {lineage: {a: freq for a, freq in annots.items() if a is not None}})

    # write data stats
    dbstat_file = os.path.join(args.outdir, f'dbstats.{timestamp}.json')

    with open(dbstat_file, 'w', encoding='utf-8') as handle:
        # total number of sequences
        NSEQS = 0
        for records in by_lineage.values():
            for variant in records.values():
                NSEQS += len(variant)  # number of samples
        val = {
            'lastupdate': timestamp.split('T')[0],
            'noseqs': NSEQS,
            'lineages': {}
        }
        for lineage, records in by_lineage.items():
            # Check if stats for the lineage are in the database
            if args.use_db and (lineage not in UPDATED_LINEAGES):
                CUR.execute(
                    f"SELECT dbstats_data FROM DBSTATS WHERE lineage = '{lineage}'")
                lineage_stats = CUR.fetchone()
                if lineage_stats is None:
                    cb.callback(
                        f"Missing DBSTATS record for lineage {lineage}",
                        level='ERROR')
                    sys.exit()
                val['lineages'][lineage] = lineage_stats['dbstats_data']
                continue

            prefix = lineage.split('.')[0]

            # resolve PANGO prefix aliases
            lname = lineage
            if (lineage.lower() not in ['unclassifiable', 'unassigned']
                    and not prefix.startswith('X')
                    and aliases[prefix] != ''):
                lname = lineage.replace(prefix, aliases[prefix])

            samples = unpack_records(records)
            ndiffs = [len(x['diffs']) for x in samples]
            val['lineages'][lineage] = {
                'nsamples': len(samples),
                'lastcoldate': max(x['covv_collection_date'] for x in samples),
                'residual': residuals[lineage] if lineage in residuals else 0,
                'max_ndiffs': max(ndiffs),
                'mean_ndiffs': sum(ndiffs) / len(ndiffs),
                'mutations': mutations[lineage],
                'infections': infection_prediction[lineage],
                'raw_lineage': lname
            }

            if args.use_db:
                CUR.execute('''
                INSERT INTO DBSTATS
                VALUES (%s, %s)
                ON CONFLICT (lineage) DO UPDATE
                SET dbstats_data = %s
                ''', [lineage, json.dumps(val['lineages'][lineage]),
                      json.dumps(val['lineages'][lineage])])

        json.dump(val, handle)

    if args.use_db:
        conn.commit()
        conn.close()

    # upload output files to webserver, requires SSH key credentials
    if not args.dry_run:
        SERVER_ROOT = 'filogeneti.ca:/var/www/html/covizu/data'
        subprocess.check_call(['scp', nwk_file, f'{SERVER_ROOT}/timetree.nwk'])
        subprocess.check_call(['scp', xbb_file, f'{SERVER_ROOT}/xbbtree.nwk'])
        subprocess.check_call(
            ['scp', clust_file, f'{SERVER_ROOT}/clusters.json'])
        subprocess.check_call(
            ['scp', dbstat_file, f'{SERVER_ROOT}/dbstats.json'])

        # upload files to EpiCoV server
        SERVER_EPICOV = 'filogeneti.ca:/var/www/html/epicov/data'
        subprocess.check_call(
            ['scp', xbb_file, f'{SERVER_EPICOV}/xbbtree.nwk'])
        subprocess.check_call(
            ['scp', nwk_file, f'{SERVER_EPICOV}/timetree.nwk'])
        subprocess.check_call(
            ['scp', dbstat_file, f'{SERVER_EPICOV}/dbstats.json'])

        # modify clusters JSON
        with open(clust_file, 'r', encoding='utf-8') as epifile:
            epicov_data = gisaid_utils.convert_json(epifile, args.infile)
            with NamedTemporaryFile('w', delete=False, encoding='utf-8') as fp:
                json.dump(epicov_data, fp=fp)  # serialize to temp file

                subprocess.check_call(
                    ['scp', fp.name, f'{SERVER_EPICOV}/clusters.json'])

    cb.callback("All done!")
