import functools
import json
import os
import re
import sys
import warnings


# define path to working directory
path = "./" # path to cloned pokay directory
text_indir = path + "data" # path to ./data folder in pokay directory
mut_dict = {}


# sub function to sort constellations by position
def constellation_comp(constellations):
    # position-based or fall back toa lexical if same position
    def compare_position(x, y):
        pos_x = re.search("^.*?(\d+)", x).group(1)
        pos_y = re.search("^.*?(\d+)", y).group(1)
        return (pos_x > pos_y) if pos_x != pos_y else (x > y)
    return sorted(constellations, key = functools.cmp_to_key(compare_position))


# sub function to match sentence
def annotations_for_sentence_match(sentence, lit_references, lit_link_evidence_summary):
    sentence_match = sentence.group(1)
    lit_link = re.search("^(.+?)\s+(https?://\S+)", sentence_match, re.MULTILINE)
    summary_text = re.search("^[A-Z](?=[a-z])", sentence_match)
    summary_text_end = re.search("\.\s*$", sentence_match)

    # lit link
    if lit_link:
        authors, link = lit_link.groups()
        lit_references.append([lit_link_evidence_summary, authors, link])
        lit_link_evidence_summary = ""
    
    # otherwise assume part of summary evidence text for lit link to come.
    # maybe a missing period between lines of evidence from same paper.
    elif summary_text:
        if not lit_link_evidence_summary: return -1
        elif re.compile("\.\s*$").match(lit_link_evidence_summary): return -1
        else: lit_link_evidence_summary += ".<br/>{}".format(summary_text)
    
    # keep the line ending
    elif summary_text_end:
        lit_link_evidence_summary += " {}<br/>".format(summary_text_end)
    
    # assume it's a sentence that runs over two lines
    else: lit_link_evidence_summary += " {}".format(sentence_match)
    
    lit_link_evidence_summary = re.sub("\.\s*<br\/>\s*\.<br\/>", ".<br\/>", sentence_match)
    return 0


# sub function to generate constellation dictionary
def generate_constellation_keys(lines, gene, qualified_txt_file, short_effect_desc, constellations_to_categories):
    lit_references = []
    lit_link_evidence_summary = ""

    for index, line in enumerate(lines):
        keys_match = re.search("^[A-Z]+\d+(?:[A-Z]+|del)", line)
        sentence = re.search("^\s*#\s*(.*)", line)
        if keys_match:
            line = line.replace(" ", "") # remove any errant whitespace
            line = line.rstrip()
            
            # sort the constellation key atoms so the naming is consistent even if the file encoding of the combos was not
            sorted_constellations = ";".join(constellation_comp(line.split(";")))
            constellation_key = "{}:{}".format(gene, sorted_constellations)

            if not lit_references:
                warnings.warn("Ignoring variant '{}' wihout literature annotations ({} line #{})".format(line, qualified_txt_file, index))
            
            if not constellation_key in constellations_to_categories:
                constellations_to_categories[constellation_key] = {}
            if not short_effect_desc in constellations_to_categories[constellation_key]:
                constellations_to_categories[constellation_key][short_effect_desc] = []
            constellations_to_categories[constellation_key][short_effect_desc].extend(lit_references)
            lit_references = []

        if sentence:
            val = annotations_for_sentence_match(sentence, lit_references, lit_link_evidence_summary)
            if val == -1: continue

        elif not re.compile("^\s*$").match(line):
            warnings.warn("Ignoring line that is not a comment starting with '#', nor a properly formatted variant ({} line #{})".format(qualified_txt_file, index))
            if lit_link_evidence_summary:
                warnings.warn("Dumping preceding comments.")
                lit_link_evidence_summary = ""
        # else it's blank
    return constellations_to_categories


# sub function to get substitution matches
def get_subs_matches(constellation, gene):
    substitutions = []
    substitution_matches = re.finditer("([A-Z])(\d+)([A-Z])", constellation)
    for sub_match in substitution_matches:
        substitution_match = re.search("([A-Z])(\d+)([A-Z])", sub_match.group())
        refAA, codon, queryAA = substitution_match.groups()
        substitutions.append({ "refAA": "{}".format(refAA), "queryAA": "{}".format(queryAA), "codon": int(codon), "gene": "{}".format(gene) })
    return substitutions
    

# sub function to get deletion matches
def get_del_matches(constellation, gene):
    deletions = []
    deletion_matches = re.finditer("([A-Z])(\d+)del", constellation)
    for del_match in deletion_matches:
        deletion_match = re.search("([A-Z])(\d+)del", del_match.group())
        refAA, codon = deletion_match.groups()
        deletions.append({ "refAA": "{}".format(refAA), "codon": int(codon), "gene": "{}".format(gene)})
    return deletions


# main function to convert text annotations from pokay data directory into JSON object
def pokay_mut_annotations():
    try:
        txt_indir = os.listdir(text_indir)
    except Exception as err:
        sys.exit("Cannot open directory {} for reading: {}".format(text_indir, err))

    constellations_to_categories = {} # { unique var combo name => { short effect desc => [lit references] } }
    for txt_file in txt_indir:
        if txt_file.startswith('.'): continue # hidden files
        
        result = re.search("^(\S+?)_(.+)\.txt$", txt_file) # bona fide text file with gene_short_effect_desc.txt naming convention
        if not result: continue
        gene, short_effect_desc = result.groups()

        qualified_txt_file = "{}/{}".format(text_indir, txt_file)
        try:
            text_file = open(qualified_txt_file)
        except Exception as err:
            sys.exit("Cannot open {} for reading: {}".format(qualified_txt_file, err))
        print("Processing {}".format(qualified_txt_file))

        lines = text_file.readlines()
        generate_constellation_keys(lines, gene, qualified_txt_file, short_effect_desc, constellations_to_categories)
        text_file.close()

    # write to JSON object
    json_entries = []
    for constellation_key in constellations_to_categories:
        gene, constellation = constellation_key.split(":", 1)
        
        substitutions = get_subs_matches(constellation, gene)
        deletions = get_del_matches(constellation, gene)        

        sorted_desc = sorted(constellations_to_categories[constellation_key])
        for short_effect_desc in sorted_desc:
            json_entries.append(
                {
                    "description": "{}".format(short_effect_desc),
                    "url": "http://people.ucalgary.ca/~gordonp/{}.html#{}".format(constellation_key, short_effect_desc),
                    "substitutions": substitutions,
                    "deletions": deletions
                }
            )
    return json_entries


# main function to reformat charaterized mutations
# new format => aa, gene, refAA, codon, queryAA
def reformatMutations(subt, dels, func):

    # define list & dictionary object
    # to store data
    mut_list = []

    # iterate through mutations
    for c, i in enumerate(subt):
        subs = i
        muts = "{}:{}:{}{}{}".format("aa",subs["gene"], subs["refAA"], subs["codon"], subs["queryAA"])
        mut_list.append(muts)

    # iterate through deletions if any
    for n, j in enumerate(dels):
        if dels == "None":
            dells = []
        else:
            dell = j
            dells = "{}:{}:{}{}{}".format("aa", dell["gene"], dell["refAA"], dell["codon"], "-")
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
outfile = open(path + 'mut_annotations.json', 'w')

# iterating through pokay generated mutation annotations
for i in range(0, len(data)):
    subs = data[i]["substitutions"]
    dels = data[i]["deletions"]
    func = data[i]["description"]

    #run function
    reformatMutations(subs, dels, func)
    
# write result to JSON file
json.dump(mut_dict, outfile, indent=4, ensure_ascii=False, separators=(',', ':'))

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