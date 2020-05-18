import argparse, os
import pandas as pd

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Performs minimap alignment to ference genome per sequence in a fasta file of genomes and collects the SNPs with respect to the reference genome to a vcf file')
    parser.add_argument('--reference_genome_path', '-r', help='path to the file holding the reference gnome',
                        required=True)
    parser.add_argument('--genomes_path', '-i',
                        help='path a fasta file holding all the covid19 genomes to be aligned to the reference',
                        required=True)
    parser.add_argument('--output_dir', '-o',
                        help='directory to hold the vcf files outputed for each pairwise alignment of a genome with the refernce genome with the collected SNPs',
                        required=True)
    parser.add_argument('--metadata_path', '-m',
                        help='path to the metadata file which can be crossed with the sequence names in the genomes path',
                        required=True)
    parser.add_argument('--start_index', '-s', help='index of sequence to start the analysis from', required=False,
                        default=0)
    parser.add_argument('--end_index', '-e', help='index of sequence to end the analysis with', required=False,
                        default=27160)

    args = parser.parse_args()
    reference_genome_path = args.reference_genome_path
    genomes_path = args.genomes_path
    output_dir = args.output_dir
    metadata_path = args.metadata_path
    start_index = int(args.start_index)
    end_index = int(args.end_index)

    concat_cmd = 'bcftools concat -o ' + output_dir + "vairants_" + str(start_index) + "_to_" + str(end_index) + ".vcf"

    # res = os.system('conda init')
    # res = os.system('conda activate CovidML')
    fasta_dir = output_dir + "genomes/"
    if not os.path.exists(fasta_dir):
        res = os.system("mkdir -p " + fasta_dir)
    vcfs_dir = output_dir + 'vcf/'
    if not os.path.exists(vcfs_dir):
        res = os.system('mkdir -p ' + vcfs_dir)
    pafs_dir = output_dir + "paf/"
    if not os.path.exists(pafs_dir):
        res = os.system("mkdir -p " + pafs_dir)

    # read the metadata file to a df
    metadata = pd.read_csv(metadata_path, sep="\t")

    # extract the sequences from the start index to the end index to distinct a fasta file, then run minimap on it, and then delete it
    with open(genomes_path, 'r') as input_file:
        lines = input_file.readlines()

    for index in range(start_index, end_index):
        line = 2 * index
        line1 = lines[line]
        sequence_name = line1.replace("\n", "").replace(">", "")
        # search for corresponding metadata for sequence name
        sequence_metadata = metadata.loc[(metadata['strain'] == sequence_name)]
        updated_sequence_name = "hCoV-19_" + sequence_name.replace("/", "_") + "_date_" + str(
            sequence_metadata['date'].values[0]).replace(" ", "_") + "_region_" + str(
            sequence_metadata['region'].values[0]).replace(" ", "_") + "_country_" + str(
            sequence_metadata['country'].values[0]).replace(" ", "_") + "_age_" + str(sequence_metadata['age'].values[0]).replace(" ", "_") + "_sex_" + str(
            sequence_metadata['sex'].values[0]).replace(" ", "_")
        line2 = lines[line + 1]
        alternative_genome_path = output_dir + '/' + updated_sequence_name + '.fa'
        with open(alternative_genome_path, 'w') as genome_file:
            genome_file.write(">" + updated_sequence_name + "\n")
            genome_file.write(line2)

        # execute minimap on the genome_path
        pas_path = pafs_dir + '/' + updated_sequence_name + '.pas'
        vcf_path = vcfs_dir + '/' + updated_sequence_name + '.vcf'
        res = os.system(
            'minimap2 -cx asm20 --cs ' + reference_genome_path + ' ' + alternative_genome_path + ' > ' + pas_path)
        res = os.system(
            'sort -k6,6 -k8,8n ' + pas_path + ' | paftools.js call -f ' + reference_genome_path + ' -L10000 -l1000 - > ' + vcf_path)

        # delete the fasta, pas and vcf files files
        res = os.system('rm -r ' + alternative_genome_path)
        res = os.system('rm -r ' + pas_path)

        # compress the variants file
        res = os.system("bgzip -c " + vcf_path + " > " + vcf_path + ".gz")
        res = os.system("bcftools index " + vcf_path + ".gz")
        res = os.system("rm -r " + vcf_path)
        vcf_path += ".gz"

        # chain th path to the vcf file to the merge command of bcf tools
        concat_cmd = concat_cmd + ' ' + vcf_path

    res = os.system(concat_cmd)
    res = os.system('rm -r ' + vcfs_dir)


