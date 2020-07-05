import argparse, sys, re, pandas as pd
from Bio import SeqIO, Entrez, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Extracts trait data according to choice of trait and a given data file, based on vlookup on genbank IDs')
    parser.add_argument('--trait_input_path', '-i', help='path to data file with trait data',
                        required=True)
    parser.add_argument('--sequence_input_path', '-i', help='path to data file with sequence genbank IDs',
                        required=True)
    parser.add_argument('--is_protein_data', '-p', help='0 if data is genomes / CDS, else 1. Default is 0', required=False, default=0)
    parser.add_argument('--trait', '-t', help='trait of choice: 0: Group, 1: Host, 2: Vector, 3: Country. Default is 0: Group', required=False, default=0)
    parser.add_argument('--output_dir', '-o', help='directory that will hold the output of the extraction: a fasta file with numerical states, a csv file that maps numerical states to strings, and a csv file that maps the genbank accessions to strings',
                        required=True)

    args = parser.parse_args()
    trait_input_path = args.trait_input_path
    sequence_input_path = args.sequence_input_path
    is_protein_data = int(args.is_protein_data)
    trait = int(args.trait)
    output_dir = args.output_dir

    # read trait input path
    trait_data = pd.read_excel(trait_input_path)

    # read records from sequence input path
    records = [record for record in SeqIO.parse(sequence_input_path, "fasta")]

    # for each record, map record name to genbank id in a dataframe that will later on hold the respective trait state
    genbank_id_regex = re.compile("gb\:(.*?)[\:|\|]", re.MULTILINE | re.DOTALL)
    nuc_id_regex = re.compile('(.*?)[\.|\:]', re.MULTILINE | re.DOTALL)
    map_df = pd.DataFrame(columns=['name', 'genbank_id', 'group', 'host', 'vector', 'country'])
    for record in records:
        name = record.description
        genbank_id = genbank_id_regex.search(name).group(1)
        # if the data is protein data, the respective protein id will not be available in trait data and must be mapped to its CDS / genome ID
        if is_protein_data:
            handle = Entrez.efetch(db="protein", id=genbank_id, rettype="gb", retmode="xml")
            record = [r for r in Entrez.parse(handle)][0]
            features = record['GBFeature_quals']
            for feature in features:
                if feature['GBQualifier_name'] == 'coded_by':
                    genbank_id = nuc_id_regex.search(feature['GBQualifier_value']).group(1)
        # get the trait data corresponding to the genbank id
        genbank_trait_data = trait_data[trait_data['GenBank Accession'].str.contains(genbank_id)]
        # get the group based on species name
        species = genbank_trait_data[['Virus Species']]