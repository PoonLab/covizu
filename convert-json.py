import argparse
import lzma
import json
import sys


def convert_json(infile, provision):
    """
    Old style contains only name, accession and coldate.
    From:  ['coldate', 'accession', 'name']
    To:    ['name', 'accession', 'location', 'coldate', 'gender', 'age', 'status']
    """
    # first generate dictionary from provision keyed by accession
    metadata = {}
    with lzma.open(provision, 'rb') as handle:
        for line in handle:
            record = json.loads(line)
            metadata.update({record['covv_accession_id']: {
                'name': record['covv_virus_name'],
                'location': record['covv_location'],
                'coldate': record['covv_collection_date'],
                'gender': record['covv_gender'],
                'age': record['covv_patient_age'],
                'status': record['covv_patient_status']
            }})

    clusters = json.load(infile)
    for cluster in clusters:
        for variant, samples in cluster['nodes'].items():
            revised = []
            for coldate, accn, name in samples:
                md = metadata.get(accn, None)
                if md is None:
                    print("Failed to retrieve metadata for accession {}".format(accn))
                    sys.exit()
                revised.append([name, accn, md['location'], coldate, md['gender'], md['age'], md['status']])

            # replace list of samples
            cluster['nodes'][variant] = revised

    return clusters


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert JSON from main branch to epicov formats.")
    parser.add_argument('json', type=argparse.FileType('r'),
                        help='input, cluster JSON file to convert')
    parser.add_argument('xz', type=str, help='input, path to xz-compressed provisioning file')
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help='output, path to write converted JSON')
    args = parser.parse_args()
    clusters = convert_json(args.json, args.xz)
    json.dump(clusters, fp=args.outfile)
