"""
script to convert JSON from main brnach to epicov format
"""
import argparse
import json
from covizu.utils.gisaid_utils import convert_json

# command line interface
parser = argparse.ArgumentParser(description="Convert JSON from main branch to epicov formats.")
parser.add_argument('json', type=argparse.FileType('r'),
                    help='input, cluster JSON file to convert')
parser.add_argument('xz', type=str, help='input, path to xz-compressed provisioning file')
parser.add_argument('outfile', type=argparse.FileType('w'),
                    help='output, path to write converted JSON')
args = parser.parse_args()

clusters = convert_json(args.json, args.xz)
json.dump(clusters, fp=args.outfile)
