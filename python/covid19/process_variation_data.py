import pandas as pd
import argparse, re, os, vcf

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(description='Process variation data from a file to produce csv of SNPs frequencies by country')
    parser.add_argument('--input_path', '-i', help='path to the variation data file in vcf format', required=True)
    parser.add_argument('--output_path', '-o', help='path to the output collecting SNP frequencies by countries in csv format', required=True)
    parser.add_argument('--metadata_to_group_by', '-m', help='metadata field to group by the SNPs frequencies by. Options: region, country, age, sex. default is country', required=False, default='country')

    args = parser.parse_args()
    input_path = args.input_path
    output_path = args.output_path
    metadata_to_group_by = args.metadata_to_group_by

    # read input file
    vcf_reader = vcf.Reader(open(input_path))
    records = [record for record in vcf_reader]

    # collect the SNPs and their frequencies across all the samples
    SNPs_to_frequency = dict()
    total_count = 0
    for record in records:
        total_count += 1
        if (str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT)) not in SNPs_to_frequency:
            SNPs_to_frequency[(str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT))] = 1
        else:
            SNPs_to_frequency[(str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT))] += 1
    for SNP in SNPs_to_frequency:
        SNPs_to_frequency[SNP] /= total_count

    # search for SNPs with at least 5% frequency across worldwide population
    frequent_SNPs = []
    for SNP in SNPs_to_frequency:
        if SNPs_to_frequency[SNP] > 0.05:
            frequent_SNPs.append(SNP)

    # now, collect a the frequencies of the filtered SNPs across countries into a dataframe
    metadata_to_group_by_regex = re.compile('_' + metadata_to_group_by + '_(.*?)[_|;]', re.DOTALL)
    data = dict()
    metadata_to_samples_count = dict()
    for record in records:
        metadata_value = metadata_to_group_by_regex.search(record.INFO['QNAME']).group(1)
        if not metadata_value in metadata_to_samples_count:
            metadata_to_samples_count[metadata_value] = 1
        else:
            metadata_to_samples_count[metadata_value] += 1
        if not metadata_value in data:
            data[metadata_value] = dict()
            for SNP in frequent_SNPs:
                data[metadata_value][SNP] = 0
        if (str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT)) in frequent_SNPs:
            data[metadata_value][(str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT))] += 1
    for metadata_value in  metadata_to_samples_count:
        for SNP in frequent_SNPs:
            data[metadata_value][SNP] /= metadata_to_samples_count[metadata_value]

    # insert data to a dataframe
    total_samples_count = sum(metadata_to_samples_count[metadata_value] for metadata_value in metadata_to_samples_count)
    df = pd.DataFrame(columns=[metadata_to_group_by, 'frequency_across_samples'] + [
        'chrom_' + str(frequent_SNP[0]) + '_pos_' + frequent_SNP[1] + '_ref_' + frequent_SNP[2] + '_alt_' + frequent_SNP[3] for frequent_SNP in frequent_SNPs])
    for metadata_value in metadata_to_samples_count:
        value = {metadata_to_group_by: metadata_value}
        value['frequency_across_samples'] = metadata_to_samples_count[metadata_value] / total_samples_count
        for SNP in frequent_SNPs:
            value['chrom_' + SNP[0] + '_pos_' + SNP[1] + '_ref_' + SNP[2] + '_alt_' + SNP[3]] = data[metadata_value][SNP]
        df = df.append(value, ignore_index = True)
    df.to_csv(output_path)

