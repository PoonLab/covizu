import os
import subprocess
from Bio import Phylo
import math
from covizu import clustering, beadplot
import sys
import json
from yaml import dump
import csv
import tempfile
import covizu.treetime
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import covizu
import math

def unpack_records(records):
    """
    by_lineage is a nested dict with the inner dicts keyed by serialized
    mutation sets (diffs).  This function is used to reconstitute the mutations
    as a list of tuples, restoring the diffs entry for each record, and to return
    a list of dicts.
    Used in:
    - treetime:retrieve_genomes()
    - get_mutations()
    - batch.py, command line interface
    See issue: https://github.com/PoonLab/covizu/issues/387

    @param records:  dict, sets of genomes under a single lineage classification,
                     each set keyed by their shared set of mutations (diffs)
    @return:  a list of dicts, each dict representing a genome and its metadata
    """
    unpacked = []
    for key, variant in records.items():
        # reconstitute the mutations defining this variant
        diffs = []
        for mutation in key.split(','):
            typ, pos, alt = mutation.split('|')
            if typ == '-':
                alt = int(alt)  # number of nucleotides in indel
            diffs.append(tuple([typ, int(pos), alt]))

        for sample in variant:
            sample.update({'diffs': diffs})
            unpacked.append(sample)
    return unpacked


def build_timetree(by_lineage, args, callback=None):
    """ Generate time-scaled tree of Pangolin lineages """

    if callback:
        callback("Parsing Pango lineage designations")
    handle = open(args.lineages)
    header = next(handle)
    if header != 'taxon,lineage\n':
        if callback:
            callback("Error: {} does not contain expected header row 'taxon,lineage'".format(args.lineages))
        sys.exit()
    lineages = {}
    for line in handle:
        try:
            taxon, lineage = line.strip().split(',')
            if taxon and lineage:
                lineages.update({taxon: lineage})
            else:
                if callback:
                    callback("Warning '{}': taxon or lineage is missing".format(line), level='WARN')
        except:
            if callback:
                callback("Warning: There is an issue with the line '{}' in lineages.csv".format(line), level='WARN')

    if callback:
        callback("Identifying lineage representative genomes")
    fasta = covizu.treetime.retrieve_genomes(by_lineage, known_seqs=lineages, ref_file=args.ref,
                                      earliest=True)

    if callback:
        callback("Reconstructing tree with {}".format(args.ft2bin))
    nwk = covizu.treetime.fasttree(fasta, binpath=args.ft2bin)

    if callback:
        callback("Reconstructing time-scaled tree with {}".format(args.ttbin))
    nexus_file = covizu.treetime.treetime(nwk, fasta, outdir=args.outdir, binpath=args.ttbin,
                                   clock=args.clock, verbosity=0)

    # writes output to treetime.nwk at `nexus_file` path
    return covizu.treetime.parse_nexus(nexus_file, fasta)


def beadplot_serial(lineage, features, args, callback=None): # pragma: no cover
    """ Compute distance matrices and reconstruct NJ trees """
    # bootstrap sampling and NJ tree reconstruction, serial mode
    trees, labels = clustering.build_trees(features, args, callback=callback)
    if trees is None:
        # lineage only has one variant, no meaningful tree
        beaddict = {'lineage': lineage, 'nodes': {}, 'edges': []}

        # use earliest sample as variant label
        intermed = [label.split('|')[::-1] for label in labels['0']]
        intermed.sort()
        variant = intermed[0][1]
        beaddict.update({'sampled_variants': len(labels)})
        beaddict['nodes'].update({variant: []})

        for coldate, accn, location, label1 in intermed:
            beaddict['nodes'][variant].append([coldate, accn, label1])
        return beaddict

    # generate majority consensus tree
    ctree = clustering.consensus(iter(trees), cutoff=args.boot_cutoff)

    # collapse polytomies and label internal nodes
    label_dict = dict([(str(idx), lst) for idx, lst in enumerate(labels)])
    ctree = beadplot.annotate_tree(ctree, label_dict, callback=callback)

    # convert to JSON format
    beaddict = beadplot.serialize_tree(ctree)
    beaddict.update({'lineage': lineage})
    beaddict.update({'sampled_variants': len(labels)})
    return beaddict


def import_labels(handle, callback=None): # pragma: no cover
    """ Load map of genome labels to tip indices from CSV file """
    result = {}
    _ = next(handle)  # skip header line
    for line in handle:
        try:
            qname, idx = line.strip('\n').split(',')
        except ValueError:
            if callback:
                callback("import_labels() failed to parse line {}".format(line), level="ERROR")
            raise  # issue #206, sequence label contains delimiter

        if idx not in result:
            result.update({idx: []})
        result[idx].append(qname)
    return result


def get_shannons(n_seqs_in_nodes):
    """Find shannon's diversity from the seuqences in the tips"""
    sample_size = sum(n_seqs_in_nodes)
    shannons = 0
    simpsons = 0
    
    for n_seq in n_seqs_in_nodes:
        proportion = (n_seq/sample_size)
        shannons += proportion*math.log(proportion)
    
    shannons = -shannons
    
    return (shannons)


def manage_collapsed_nodes(labels, tree):
    """Add collapsed node keys to the labels dictionary"""
    new_labels = labels
    for clade in tree.get_terminals() + tree.get_nonterminals():
        name = clade.name
        if name == None:
            continue
        elif '|' in name:
            combined_list = []
            for title in name.split('|'):
                combined_list = combined_list + labels[title]
            new_labels[name] = combined_list
    return new_labels


def get_tree_summary_stats(tree, Ne, label_dict, keep_temp = False):
    """Write out file of summary stats including number of unsampled lineages,
    diversity metrics and root to tip regression"""
    internal_nodes = tree.get_nonterminals()
    terminal_nodes = tree.get_terminals()
   
    #Unsampled nodes are nodes with no sequence ID associated with them
    num_terminal_nodes = len(terminal_nodes)
    seqs_in_term_nodes = [len(label_dict[node.name]) for node in terminal_nodes]
    seqs_in_internal_nodes = []
    for node in internal_nodes:
        if node.name != None:
            seqs_in_internal_nodes = seqs_in_internal_nodes + [len(label_dict[node.name])]
        
    #Calculate unsampled lineages and shannon's diversity
    unsampled_count = sum([node.name==None for node in internal_nodes])
    
    diversity = get_shannons(seqs_in_term_nodes)
    
    #Create dictionary of summary stats
    summary_stats = {'unsampled_lineage_count': unsampled_count,
                         'shannons_diversity': diversity,
                         'Ne': Ne
                         }
    return summary_stats


def find_Ne(tree, labels_filename, sequence_labels, keep_temp = False):
    """Run beta skyline estimation implemented into R"""
    
    if (keep_temp):
        tree_filename = "combined_tree.tree"
    else:
        #Make a temporary file containing the tree
        temp_tree, tree_filename = tempfile.mkstemp(suffix =".tree")

    Phylo.write(tree, tree_filename, "nexus")
    
    #Get the value of Ne
    Ne = subprocess.Popen("Rscript " + os.path.dirname(__file__) + "/skyline_est.R -t " + 
            tree_filename + " -l " + labels_filename, shell = True, 
            stdout=subprocess.PIPE).stdout.read().decode() 
    

    if (not keep_temp):
        #Remove the temporary file
        os.remove(tree_filename)
    
    return Ne


def get_diversity(indexed, labels):
    """
    Calculate an analogue to the nucleotide diversity (the expected number of
    differences between two randomly sampled genomes).
    :param indexed:  list, sets of feature indices for each variant
    :param labels:  dict, {variant number: [sequence names]}
    """
    nvar = len(indexed)
    counts = [len(v) for v in labels]  # number of genomes per variant
    total = sum(counts)
    result = 0
    for i in range(0, nvar-1):
        fi = counts[i] / total  # frequency of i-th variant
        for j in range(i+1, nvar):
            fj = counts[j] / total
            ndiff = len(indexed[i] ^ indexed[j])  # symmetric difference
            result += 2*ndiff * fi * fj
    return (result * (total / (total-1)))


def make_beadplots(by_lineage, args, callback=None, t0=None, txtfile='minor_lineages.txt',
                   recode_file="recoded.json"):
    """
    Wrapper for beadplot_serial - divert to clustering.py in MPI mode if
    lineage has too many genomes.

    :param by_lineage:  dict, feature vectors stratified by lineage
    :param args:  Namespace, from argparse.ArgumentParser()
    :param t0:  float, datetime.timestamp.
    :param txtfile:  str, path to file to write minor lineage names
    :param recode_file:  str, path to JSON file to write recoded lineage data

    :return:  list, beadplot data by lineage
    """

    # recode data into variants and serialize
    if callback:
        callback("Recoding features, compressing variants..")
    recoded = {}
    for lineage, records in by_lineage.items():
        union, labels, indexed = clustering.recode_features(records, limit=args.max_variants)

        # serialize tuple keys (features of union), #335
        union = dict([("{0}|{1}|{2}".format(*feat), idx) for feat, idx in union.items()])
        indexed = [list(s) for s in indexed]  # sets cannot be serialized to JSON, #335
        recoded.update({lineage: {'union': union, 'labels': labels,
                                  'indexed': indexed}})

    with open(recode_file, 'w') as handle:
        json.dump(recoded, handle)

    # partition lineages into major and minor categories
    intermed = [(len(features), lineage) for lineage, features in by_lineage.items()
                if len(features) < args.mincount]
    intermed.sort(reverse=True)  # descending order
    minor = dict([(lineage, None) for _, lineage in intermed if lineage is not None])

    # export minor lineages to text file
    with open(txtfile, 'w') as handle:
        for lineage in minor:
            handle.write('{}\n'.format(lineage))

    # launch MPI job across minor lineages
    if callback:
        callback("start MPI on minor lineages")
    cmd = ["mpirun", "--machinefile", args.machine_file, "python3", "covizu/clustering.py",
           recode_file, txtfile,  # positional arguments <JSON file>, <str>
           "--mode", "flat",
           "--max-variants", str(args.max_variants),
           "--nboot", str(args.nboot),
           "--outdir", args.outdir,
           "--binpath", args.binpath  # RapidNJ
           ]
    if t0:
        cmd.extend(["--timestamp", str(t0)])
    subprocess.check_call(cmd)

    # process major lineages
    for lineage, features in by_lineage.items():
        if lineage in minor:
            continue

        if callback:
            callback('start {}, {} entries'.format(lineage, len(features)))

        cmd = [
            "mpirun", "--machinefile", args.machine_file, "python3", "covizu/clustering.py",
            recode_file, lineage,  # positional arguments <JSON file>, <str>
            "--mode", "deep",
            "--max-variants", str(args.max_variants),
            "--nboot", str(args.nboot),
            "--outdir", args.outdir,
            "--binpath", args.binpath
        ]
        if t0:
            cmd.extend(["--timestamp", str(t0)])
        subprocess.check_call(cmd)

    # parse output files
    if callback:
        callback("Parsing output files")

    result = []
    inf_predict = {}

    # Load required R packages
    tidyquant = importr('tidyquant')
    yaml = importr('yaml')
    matrixStats = importr('matrixStats')

    # Read Models
    robjects.r('increasing_mods <- readRDS("{}")'.format(os.path.join(covizu.__path__[0], "hunepi/infections_increasing_model_comparisons.rds")))
    robjects.r('infections_mods <- readRDS("{}")'.format(os.path.join(covizu.__path__[0], "hunepi/num_infections_model_comparisons.rds")))

    # Function to make estimates from each model
    robjects.r('''
    estimate_vals <- function(models, predict_dat, exp = FALSE){
        prediction_df <- data.frame(sapply(models, predict, newdata = predict_dat, type = "response"))
        if (exp) {
            prediction_df <- exp(prediction_df)
        }
        return(prediction_df)
    }
    ''')

    for lineage in recoded:
        # import trees
        lineage_name = lineage.replace('/', '_')  # issue #297
        outfile = open('data/{}.nwk'.format(lineage_name))
        trees = Phylo.parse(outfile, 'newick')  # note this returns a generator

        label_dict = recoded[lineage]['labels']

        if len(label_dict) == 1:
            # handle case of only one variant
            # lineage only has one variant, no meaningful tree
            beaddict = {'nodes': {}, 'edges': []}

            # use earliest sample as variant label
            intermed = [label.split('|')[::-1] for label in label_dict['0']]
            intermed.sort()
            variant = intermed[0][1]
            beaddict['nodes'].update({variant: []})

            for coldate, accn, location, label1 in intermed:
                beaddict['nodes'][variant].append([coldate, accn, location, label1])
            
            inf_predict.update({lineage: 0})
        else:
            # generate beadplot data
            ctree = clustering.consensus(trees, cutoff=args.boot_cutoff, callback=callback)
            outfile.close()  # done with Phylo.parse generator

            # incorporate hunipie
            clabel_dict = manage_collapsed_nodes(label_dict, ctree)

            # Write labels to file for Ne estimation
            with open('clabels_file.csv', 'w') as csv_file:
                writer = csv.writer(csv_file)
                for key, value in clabel_dict.items():
                    writer.writerow([key, value])

            cne = find_Ne(ctree, 'clabels_file.csv', False)

            # Collapse tree and manage the collapsed nodes
            tree = beadplot.collapse_polytomies(ctree)
            clabel_dict = manage_collapsed_nodes(label_dict, tree)
            summary_stats = get_tree_summary_stats(tree, cne, clabel_dict, False)

            indexed = [set(l) for l in recoded[lineage]['indexed']]
            pi = get_diversity(indexed, label_dict)
            summary_stats['pi'] = pi
            summary_stats['sample_size'] = len(by_lineage[lineage])

            if cne == '':
                summary_stats['Ne'] = 'NaN'

            with open('{}.yml'.format(lineage), 'w') as yml:
                #Write out files from the script
                dump(summary_stats, yml)
          
            # Read in yaml data
            robjects.r('sum_stat_dat <- as.data.frame(read_yaml("{}.yml"))'.format(lineage))
            robjects.r('sum_stat_dat$Ne <- as.numeric(sum_stat_dat$Ne)')

            robjects.r('''
            increasing_predict_prob <- estimate_vals(increasing_mods, sum_stat_dat)
            
            if(!is.nan(sum_stat_dat$Ne)) {
                pred_prob <- increasing_predict_prob[which(rownames(increasing_predict_prob) == "HUNePi.1"), ]
            } else {
                pred_prob <- increasing_predict_prob[which(rownames(increasing_predict_prob) == "HUPi.1"), ]
            }

            predicted_increase <- pred_prob > 0.5

            if(predicted_increase){
                infections_predictions <- 
                    estimate_vals(infections_mods, sum_stat_dat, exp = T)
                if(!is.nan(sum_stat_dat$Ne)) {
                    predicted_infections <- infections_predictions[which(rownames(infections_predictions) == "HUNePi.1"), ]
                } else {
                    predicted_infections <- infections_predictions[which(rownames(infections_predictions) == "HUPi.1"), ]
                }
            } else {
                predicted_infections <- -1
            }
            ''')

            predicted_infections = list(robjects.r('predicted_infections'))[0]
            inf_predict.update({lineage: predicted_infections})

            # csummary_stats = get_tree_summary_stats(ctree, cne, clabel_dict, False)

            ctree = beadplot.annotate_tree(ctree, label_dict)
            beaddict = beadplot.serialize_tree(ctree)

            # beaddict = beaddict | csummary_stats
        
        outfile.close()  # done with Phylo.parse generator
        beaddict.update({'sampled_variants': len(label_dict)})
        beaddict.update({'lineage': lineage})
        result.append(beaddict)

    return result, inf_predict


def get_mutations(by_lineage):
    """
    Extract common mutations from feature vectors for each lineage
    :param by_lineage:  dict, return value from process_feed()
    :return:  dict, common mutations by lineage
    """
    result = {}
    for lineage, records in by_lineage.items():
        samples = unpack_records(records)

        # enumerate features
        counts = {}
        for sample in samples:
            for diff in sample['diffs']:
                feat = tuple(diff)
                if feat not in counts:
                    counts.update({feat: 0})
                counts[feat] += 1

        # filter for mutations that occur in at least half of samples
        common = dict([(feat, count/len(samples)) for feat, count in counts.items()
                       if count/len(samples) >= 0.5])
        result.update({lineage: common})

    return result

