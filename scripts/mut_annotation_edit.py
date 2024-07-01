"""
what mut_annotation_edit.py does
"""
import functools
import json
import logging
import os
import re
import sys


# define path to working directory
PATH = "./"  # path to cloned pokay directory
TEXT_INDIR = PATH + "data"  # path to ./data folder in pokay directory
mut_dict = {}

def constellation_comp(constellations):
    """sub function to sort constellations by position"""
    # position-based or fall back toa lexical if same position
    def compare_position(x_coord, y_coord):
        pos_x = re.search("^.*?(\\d+)", x_coord).group(1)
        pos_y = re.search("^.*?(\\d+)", y_coord).group(1)
        return (pos_x > pos_y) if pos_x != pos_y else (x_coord > y_coord)
    return sorted(constellations, key=functools.cmp_to_key(compare_position))

def annotations_for_sentence_match(
        sentence,
        lit_references,
        lit_link_evidence_summary):
    """sub function to match sentence"""
    sentence_match = sentence.group(1)
    lit_link = re.search(
        "^(.+?)\\s+(https?://\\S+)",
        sentence_match,
        re.MULTILINE)
    summary_text = re.search("^[A-Z](?=[a-z])", sentence_match)
    summary_text_end = re.search("\\.\\s*$", sentence_match)

    # lit link
    if lit_link:
        authors, link = lit_link.groups()
        lit_references.append([lit_link_evidence_summary, authors, link])
        lit_link_evidence_summary = ""

    # otherwise assume part of summary evidence text for lit link to come.
    # maybe a missing period between lines of evidence from same paper.
    elif summary_text:
        if not lit_link_evidence_summary:
            return -1
        if re.compile("\\.\\s*$").match(lit_link_evidence_summary):
            return -1
        lit_link_evidence_summary += f".<br/>{summary_text}"

    # keep the line ending
    elif summary_text_end:
        lit_link_evidence_summary += f" {summary_text_end}<br/>"

    # assume it's a sentence that runs over two lines
    else:
        lit_link_evidence_summary += f" {sentence_match}"

    lit_link_evidence_summary = re.sub(
        "\\.\\s*<br\\/>\\s*\\.<br\\/>",
        ".<br\\/>",
        sentence_match)
    return 0

def generate_constellation_keys(
        lines,
        gene,
        qualified_txt_file,
        short_effect_desc,
        constellations_to_categories):
    """sub function to generate constellation dictionary"""
    lit_references = []
    lit_link_evidence_summary = ""

    for index, line in enumerate(lines):
        keys_match = re.search("^[A-Z]+\\d+(?:[A-Z]+|del)", line)
        sentence = re.search("^\\s*#\\s*(.*)", line)
        if keys_match:
            line = line.replace(" ", "")  # remove any errant whitespace
            line = line.rstrip()

            # sort the constellation key atoms so the naming is consistent even
            # if the file encoding of the combos was not
            sorted_constellations = ";".join(
                constellation_comp(line.split(";")))
            constellation_key = f"{gene}:{sorted_constellations}"

            if not lit_references:
                logging.warning(
                    "Ignoring variant '%s' without literature annotations (%s line #%s)",
                    line, qualified_txt_file, index)
            if constellation_key not in constellations_to_categories:
                constellations_to_categories[constellation_key] = {}
            if short_effect_desc not in constellations_to_categories[constellation_key]:
                constellations_to_categories[constellation_key][short_effect_desc] = [
                ]
            constellations_to_categories[constellation_key][short_effect_desc].extend(
                lit_references)
            lit_references = []

        elif sentence:
            val = annotations_for_sentence_match(
                sentence, lit_references, lit_link_evidence_summary)
            if val == -1:
                continue

        elif not re.compile("^\\s*$").match(line):
            logging.warning(
                "Ignoring line that is not a comment starting with '#'," 
                "nor a properly formatted variant (%s line #%s)",
                qualified_txt_file, index)
            if lit_link_evidence_summary:
                logging.warning("Dumping preceding comments.")
                lit_link_evidence_summary = ""
        # else it's blank
    return constellations_to_categories

def get_subs_matches(constellation, gene):
    """sub function to get substitution matches"""
    substitutions = []
    substitution_matches = re.finditer("([A-Z])(\\d+)([A-Z])", constellation)
    for sub_match in substitution_matches:
        substitution_match = re.search(
            "([A-Z])(\\d+)([A-Z])", sub_match.group())
        ref_aa, codon, query_aa = substitution_match.groups()
        substitutions.append({"refAA": f"{ref_aa}",
                              "queryAA": f"{query_aa}",
                              "codon": int(codon),
                              "gene": f"{gene}"})
    return substitutions

def get_del_matches(constellation, gene):
    """sub function to get deletion matches"""
    deletions = []
    deletion_matches = re.finditer("([A-Z])(\\d+)del", constellation)
    for del_match in deletion_matches:
        deletion_match = re.search("([A-Z])(\\d+)del", del_match.group())
        ref_aa, codon = deletion_match.groups()
        deletions.append(
            {"refAA": f"{ref_aa}", "codon": int(codon), "gene": f"{gene}"})
    return deletions

def pokay_mut_annotations():
    """main function to convert text annotations from pokay data directory into JSON object"""
    try:
        txt_indir = os.listdir(TEXT_INDIR)
    except FileNotFoundError as err:
        sys.exit(f"Cannot open directory {TEXT_INDIR} for reading: {err}")

    # { unique var combo name => { short effect desc => [lit references] } }
    constellations_to_categories = {}
    for txt_file in txt_indir:
        if txt_file.startswith('.'):
            continue  # hidden files

        # bona fide text file with gene_short_effect_desc.txt naming convention
        result = re.search("^(\\S+?)_(.+)\\.txt$", txt_file)
        if not result:
            continue
        gene, short_effect_desc = result.groups()

        qualified_txt_file = f"{TEXT_INDIR}/{txt_file}"
        try:
            with open(qualified_txt_file, 'r', encoding='utf-8') as text_file:
                lines = text_file.readlines()
                generate_constellation_keys(
                    lines,
                    gene,
                    qualified_txt_file,
                    short_effect_desc,
                    constellations_to_categories)
                text_file.close()
        except FileNotFoundError as err:
            sys.exit(f"Cannot open {qualified_txt_file} for reading: {err}")
        print(f"Processing {qualified_txt_file}")

    # write to JSON object
    json_entries = write_pokay_json(constellations_to_categories)
    return json_entries

def write_pokay_json(const_cat_dict):
    """take constellation to cat dict and organize it to write out to json"""
    json_ent_out = []
    for constellation_key, value in const_cat_dict:
        gene, constellation = constellation_key.split(":", 1).items()

        substitutions = get_subs_matches(constellation, gene)
        deletions = get_del_matches(constellation, gene)

        sorted_desc = sorted(value)
        for short_effect_desc in sorted_desc:
            json_ent_out.append(
                {
                    "description": f"{short_effect_desc}",
                    "url": f"http://people.ucalgary.ca/" \
                           f"~gordonp/{constellation_key}.html#{short_effect_desc}",
                    "substitutions": substitutions,
                    "deletions": deletions})
    return json_ent_out

def reformat_mutations(subt, dels, func):
    """main function to reformat charaterized mutations 
    new format => aa, gene, refAA, codon, queryAA"""
    # define list & dictionary object
    # to store data
    mut_list = []

    # iterate through mutations
    for i in subt:
        substitution = i
        muts = f"aa:{substitution['gene']}:{substitution['refAA']}" \
               f"{substitution['codon']}{substitution['queryAA']}"
        mut_list.append(muts)

    # iterate through deletions if any
    for j in dels:
        if dels == "None":
            dells = []
        else:
            dell = j
            dells = f"aa:{dell['gene']}:{dell['refAA']}{dell['codon']}-"
            mut_list.append(dells)

    # keyed by mutations
    for mut in mut_list:
        if mut not in mut_dict:
            mut_dict[mut] = [func]
        elif func not in mut_dict[mut]:
            mut_dict[mut].append(func)


# get JSON object from pokay directory
data = pokay_mut_annotations()

# output JSON
with open(PATH + 'mut_annotations.json', 'w', encoding='utf-8') as outfile:

    # iterating through pokay generated mutation annotations
    for ind, info in enumerate(data):
        out_subs = info["substitutions"]
        out_dels = info["deletions"]
        out_func = info["description"]

        # run function
        reformat_mutations(out_subs, out_dels, out_func)

    # write result to JSON file
    json.dump(
        mut_dict,
        outfile,
        indent=4,
        ensure_ascii=False,
        separators=(
            ',',
            ':'))

    # closing files
    outfile.close()


# fix JSON file issues, add commas
# and a square brackets at the begin and end of file
# with open(path + "mut_annotations.json", "r+") as fyl:
# 	edit_f = fyl.read()
# 	fyl.seek(0)
# 	j_son = '[' + edit_f  + ']'
# 	j_son = j_son.replace('}{', '},{')
# 	fyl.write(j_son)
# 	fyl.close
