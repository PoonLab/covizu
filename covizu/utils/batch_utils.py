"""batch utils docstring"""
import os
import subprocess
import math
import sys
import json
import csv
from tempfile import NamedTemporaryFile
from datetime import datetime
from Bio import Phylo
from covizu import clustering, beadplot, treetime
from rpy2 import robjects
from rpy2.robjects.packages import importr
import covizu


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
            try:
                typ, pos, alt = mutation.split('|')
            except ValueError:
                # Handle case when there are no diffs
                continue
            if typ == '-':
                alt = int(alt)  # number of nucleotides in indel
            diffs.append(tuple([typ, int(pos), alt]))

        for sample in variant:
            sample.update({'diffs': diffs})
            unpacked.append(sample)
    return unpacked

def parse_lineages(in_args, in_callback):
    """parse lineages from build_timetree"""
    with open(in_args.lineages, encoding='utf-8') as handle:
        header = next(handle)
        if header != 'taxon,lineage\n':
            if in_callback:
                in_callback(
                    f"Error: {in_args.lineages} "
                    f"does not contain expected header row 'taxon,lineage'")
            sys.exit()
        lineages = {}
        for line in handle:
            try:
                taxon, lineage = line.strip().split(',')
                if taxon and lineage:
                    lineages.update({taxon: lineage})
                else:
                    if in_callback:
                        in_callback(
                            f"Warning '{line}': taxon or lineage is missing",
                            level='WARN')
            except ValueError:
                if in_callback:
                    in_callback(
                        f"Warning: There is an issue with the line '{line}' in lineages.csv",
                        level='WARN')
        return lineages

def build_timetree(
        by_lineage,
        args,
        outgroup=None,
        callback=None,
        debug=False):
    """
    Generate time-scaled tree of Pangolin lineages
    @param by_lineage:  dict, lists of records keyed by feature sets, keyed by lineage
    @param args:  argparse.Namespace, arguments from CLI
    @param outgroup:  str, optionally specify path to FASTA file containing outgroup sequence
    @param debug:  bool, if True, write out FASTA of sequences from retrieve_genomes
    @return str, Newick tree string produced by parse_nexus
    """
    if callback:
        callback("Parsing Pango lineage designations")
    lineages = parse_lineages(args, callback)

    if callback:
        callback("Identifying lineage representative genomes")
    fasta = covizu.treetime.retrieve_genomes(
        by_lineage,
        known_seqs=lineages,
        ref_file=args.ref,
        outgroup=outgroup,
        earliest=True)
    if debug:
        with open(f"build_timetree_dump.{datetime.now().isoformat()}.fasta", 'w',
                  encoding='utf-8') as outfile:
            for header, sequence in fasta.items():
                outfile.write(f">{header}\n{sequence}\n")

    if callback:
        callback(f"Reconstructing tree with {args.ft2bin}")
    nwk = covizu.treetime.fasttree(fasta, binpath=args.ft2bin)

    if callback:
        callback(f"Reconstructing time-scaled tree with {args.ttbin}")
    nexus_file = covizu.treetime.treetime(
        nwk,
        fasta,
        outdir=args.outdir,
        binpath=args.ttbin,
        clock=args.clock,
        verbosity=0)

        # writes output to treetime.nwk at `nexus_file` path
    return covizu.treetime.parse_nexus(nexus_file, fasta, callback)


def beadplot_serial(lineage, features, args, callback=None):  # pragma: no cover
    """ Compute distance matrices and reconstruct NJ trees """
    # bootstrap sampling and NJ tree reconstruction, serial mode
    trees, labels = clustering.build_trees(features, args, callback=callback)
    if trees is None:
        # lineage only has one variant, no meaningful tree
        beaddict = {'lineage': lineage, 'nodes': {}, 'edges': []}

        # use earliest sample as variant label
        intermed = sorted([label.split('|')[::-1] for label in labels['0']])
        variant = intermed[0][1]
        beaddict.update({'sampled_variants': len(labels)})
        beaddict['nodes'].update({variant: []})

        for coldate, accn, _, label1 in intermed:
            beaddict['nodes'][variant].append([coldate, accn, label1])
        return beaddict

    # generate majority consensus tree
    ctree = clustering.consensus(iter(trees), cutoff=args.boot_cutoff)

    # collapse polytomies and label internal nodes
    label_dict = {str(idx):lst for idx, lst in enumerate(labels)}
    ctree = beadplot.annotate_tree(ctree, label_dict, callback=callback)

    # convert to JSON format
    beaddict = beadplot.serialize_tree(ctree)
    beaddict.update({'lineage': lineage})
    beaddict.update({'sampled_variants': len(labels)})
    return beaddict


def import_labels(handle, callback=None):  # pragma: no cover
    """ Load map of genome labels to tip indices from CSV file """
    result = {}
    _ = next(handle)  # skip header line
    for line in handle:
        try:
            qname, idx = line.strip('\n').split(',')
        except ValueError:
            if callback:
                callback(
                    f"import_labels() failed to parse line {line}",
                    level="ERROR")
            raise  # issue #206, sequence label contains delimiter

        if idx not in result:
            result.update({idx: []})
        result[idx].append(qname)
    return result


def get_shannons(n_seqs_in_nodes):
    """Find shannon's diversity from the seuqences in the tips"""
    sample_size = sum(n_seqs_in_nodes)
    shannons = 0

    for n_seq in n_seqs_in_nodes:
        proportion = n_seq / sample_size
        shannons += proportion * math.log(proportion)

    shannons = -shannons

    return shannons


def manage_collapsed_nodes(labels, tree):
    """Add collapsed node keys to the labels dictionary"""
    new_labels = labels
    for clade in tree.get_terminals() + tree.get_nonterminals():
        name = clade.name
        if name is None:
            continue
        if '|' in name:
            combined_list = []
            for title in name.split('|'):
                combined_list = combined_list + labels[title]
            new_labels[name] = combined_list
    return new_labels


def get_tree_summary_stats(tree, in_ne, label_dict):
    """Write out file of summary stats including number of unsampled lineages,
    diversity metrics and root to tip regression"""
    internal_nodes = tree.get_nonterminals()
    terminal_nodes = tree.get_terminals()

    # Unsampled nodes are nodes with no sequence ID associated with them

    seqs_in_term_nodes = [len(label_dict[node.name])
                          for node in terminal_nodes]
    seqs_in_internal_nodes = []
    for node in internal_nodes:
        if node.name is not None:
            seqs_in_internal_nodes = seqs_in_internal_nodes + \
                [len(label_dict[node.name])]

    # Calculate unsampled lineages and shannon's diversity
    unsampled_count = sum(node.name is None for node in internal_nodes)

    diversity = get_shannons(seqs_in_term_nodes)

    # Create dictionary of summary stats
    summary_stats = {'unsampled_lineage_count': unsampled_count,
                     'shannons_diversity': diversity,
                     'Ne': in_ne
                     }
    return summary_stats


def find_ne(tree, labels_filename):
    """Run beta skyline estimation implemented into R"""

    # Make a temporary file containing the tree
    with NamedTemporaryFile('w', delete=False) as tree_filename:
        Phylo.write(tree, tree_filename.name, "nexus")


    # Load required R packages
    ape = importr('ape')
    phytools = importr('phytools')
    LambdaSkyline = importr('LambdaSkyline')

    robjects.r.assign("tree_filename", tree_filename.name)
    robjects.r.assign("sequence_labels_file", labels_filename)

    try:
        robjects.r('''
        set.seed(123456)

        tree = read.nexus(tree_filename)
        sequence_labels = read.csv(sequence_labels_file)
        colnames(sequence_labels) = c("index", "value")

        #Adjust tree to include branches of length 0 on identical sequences
        tip_count = table(sequence_labels$index)
        add_tip_count = data.frame(tip_count - 1)

        for (tip_place_in_table in 1:nrow(add_tip_count)){
            tip_name = add_tip_count[tip_place_in_table,1]
            freq = add_tip_count[tip_place_in_table,2]
            if(freq != 0){
                for (counter in 1:freq){
                tree <- bind.tip(tree, paste0(tip_name,"_", counter), edge.length = 0,
                                    where=which(tree$tip.label == tip_name))
                }
            }
        }

        #Run skyline estimation
        alpha = betacoal.maxlik(tree)
        skyline = (skyline.multi.phylo(tree, alpha$p1))

        #Output skyline estimation
        pop_sizes <- head(skyline$population.size, n = 5)
        mean_pop_size <- mean(pop_sizes, na.rm = TRUE)
        ''')
        ne_out = list(robjects.r('mean_pop_size'))[0]
    except Exception as error:
        print(f'Error in finding Ne, {error}')
        ne_out = ''

    # Remove the temporary file
    os.remove(tree_filename.name)

    return ne_out


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
    for i in range(0, nvar - 1):
        freq_i = counts[i] / total  # frequency of i-th variant
        for j in range(i + 1, nvar):
            freq_j = counts[j] / total
            ndiff = len(indexed[i] ^ indexed[j])  # symmetric difference
            result += 2 * ndiff * freq_i * freq_j
    return (result * (total / (total - 1)))


def parse_alias(alias_file):
    """
    Parse PANGO alias_key.json file contents, excluding entries with empty string values.
    :param alias_file:  str, path to JSON file
    """
    alias = {}
    with open(alias_file, 'r', encoding='utf-8') as handle:
        alias = json.loads(handle.read())
        for key, value in alias.items():
            if value != '':
                alias.update({key: value})
    return alias

def make_beadplots(
        by_lineage,
        args,
        callback=None,
        initial_time=None,
        updated_lineages=None,
        txtfile='minor_lineages.txt',
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
        if updated_lineages is not None and lineage not in updated_lineages:
            continue
        union, labels, indexed = clustering.recode_features(
            records, limit=args.max_variants)

        # serialize tuple keys (features of union), #335
        union = {f"{feat[0]}|{feat[1]}|{feat[2]}": idx for feat, idx in union.items()}

        # sets cannot be serialized to JSON, #335
        indexed = [list(s) for s in indexed]
        recoded.update({lineage: {'union': union, 'labels': labels,
                                  'indexed': indexed}})


    with open(recode_file, 'w', encoding='utf-8') as handle:
        json.dump(recoded, handle)

    # partition lineages into major and minor categories
    intermed = [(len(features), lineage) for lineage, features in by_lineage.items()
                if len(features) < args.mincount and
                (updated_lineages is None or lineage in updated_lineages)]
    intermed.sort(reverse=True)  # descending order
    minor = {lineage: None for _, lineage in intermed if lineage is not None}



    # export minor lineages to text file
    with open(txtfile, 'w', encoding='utf-8') as handle:
        for lineage in minor:
            handle.write(f'{lineage}\n')

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
    if initial_time:
        cmd.extend(["--timestamp", str(initial_time)])
    subprocess.check_call(cmd)

    # process major lineages
    for lineage, features in by_lineage.items():
        if lineage in minor or (
                updated_lineages is not None and lineage not in updated_lineages):
            continue

        if callback:
            callback(f'start {lineage}, {len(features)} entries')

        cmd = [
            "mpirun", "--machinefile", args.machine_file, "python3", "covizu/clustering.py",
            recode_file, lineage,  # positional arguments <JSON file>, <str>
            "--mode", "deep",
            "--max-variants", str(args.max_variants),
            "--nboot", str(args.nboot),
            "--outdir", args.outdir,
            "--binpath", args.binpath
        ]
        if initial_time:
            cmd.extend(["--timestamp", str(initial_time)])
        subprocess.check_call(cmd)

    # parse output files
    if callback:
        callback("Parsing output files")

    result = []
    inf_predict = {}


    # Load required R packages
    tidyquant = importr('tidyquant')
    matrixStats = importr('matrixStats')

    path_1 = os.path.join(covizu.__path__[0], 'hunepi/infections_increasing_model_comparisons.rds')
    path_2 = os.path.join(covizu.__path__[0], 'hunepi/num_infections_model_comparisons.rds')
    # Read Models
    robjects.r(
        f"increasing_mods <- readRDS(\"{path_1}\")")
    robjects.r(
        f"infections_mods <- readRDS(\"{path_2}\")")



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

    for lineage, value in recoded.items():
        # import trees
        lineage_name = lineage.replace('/', '_')  # issue #297
        with open(f'{args.outdir}/{lineage_name}.nwk', encoding='utf-8') as outfile:
            trees = Phylo.parse(outfile, 'newick')  # returns a generator

            label_dict = recoded[lineage]['labels']

            if len(label_dict) == 1:
                # handle case of only one variant
                # lineage only has one variant, no meaningful tree
                beaddict = {'nodes': {}, 'edges': []}

                # use earliest sample as variant label
                intermed = sorted([label.split('|')[::-1]
                                for label in label_dict['0']])
                variant = intermed[0][1]
                beaddict['nodes'].update({variant: []})

                for coldate, accn, location, label1 in intermed:
                    beaddict['nodes'][variant].append(
                        [coldate, accn, location, label1])

                inf_predict.update({lineage: 0})
            else:
                # generate beadplot data
                ctree = clustering.consensus(
                    trees, cutoff=args.boot_cutoff, callback=callback)
                outfile.close()  # done with Phylo.parse generator

                # incorporate hunipie
                clabel_dict = manage_collapsed_nodes(label_dict, ctree)

                with NamedTemporaryFile('w', delete=False) as labels_filename:

                    # Write labels to file for Ne estimation
                    writer = csv.writer(labels_filename)
                    for key, value in clabel_dict.items():
                        writer.writerow([key, value])

                cne = find_ne(ctree, labels_filename.name)

                # Remove temporary file with labels
                os.remove(labels_filename.name)

                # Collapse tree and manage the collapsed nodes
                tree = beadplot.collapse_polytomies(ctree)
                clabel_dict = manage_collapsed_nodes(label_dict, tree)
                summary_stats = get_tree_summary_stats(tree, cne, clabel_dict)

                indexed = [set(l) for l in recoded[lineage]['indexed']]
                p_i = get_diversity(indexed, label_dict)
                summary_stats['pi'] = p_i
                summary_stats['sample_size'] = len(by_lineage[lineage])

                if cne == '':
                    summary_stats['Ne'] = 'NaN'

                sum_stat_dat = robjects.vectors.DataFrame(summary_stats)
                robjects.r.assign('sum_stat_dat', sum_stat_dat)

                robjects.r('''
                sum_stat_dat$Ne <- as.numeric(sum_stat_dat$Ne)
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

                ctree = beadplot.annotate_tree(ctree, label_dict)
                beaddict = beadplot.serialize_tree(ctree)

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
        common = {feat: count / len(samples) for feat,
                  count in counts.items() if count / len(samples) >= 0.5}
        result.update({lineage: common})

    return result
