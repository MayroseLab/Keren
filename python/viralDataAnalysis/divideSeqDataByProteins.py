from Bio import SeqIO
import re, os,argparse

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Divides sequences in a file by protein name')
    parser.add_argument('--input_path', '-i', help='path to the sequence data that holds sequences of different proteins / coding sequences in a fasta format',
                        required=True)
    parser.add_argument('--output_dir', '-o', help='directory that will hold the csv output of the analysis',
                        required=True)

    args = parser.parse_args()
    input_path = args.input_path
    output_dir = args.output_dir

    records = [record for record in SeqIO.parse(input_path, "fasta")]

    # get unique protein names
    protein_name_regex = re.compile("Protein Name\:(.*?)\|", re.MULTILINE | re.DOTALL)
    protein_name_to_records = dict()
    for record in records:
        try:
            protein_name = protein_name_regex.search(record.description).group(1).replace("/","_")
            if 'polyprotein' in protein_name and protein_name != 'truncated polyprotein NS1 prime':
                if protein_name in protein_name_to_records:
                    protein_name_to_records[protein_name].append(record)
                else:
                    protein_name_to_records[protein_name] = [record]
        except Exception as e:
            pass

    # for each protein name, create a sequence data file in the output directory with the respective records
    output_path = output_dir + "unaligned_polyprotein.fasta"
    records = []
    for protein_name in protein_name_to_records:
        records += protein_name_to_records[protein_name]
    SeqIO.write(records, output_path, "fasta")
    # for protein_name in protein_name_to_records:
    #     output_path = output_dir + protein_name + ".fasta"
    #     SeqIO.write(protein_name_to_records[protein_name], output_path, "fasta")


