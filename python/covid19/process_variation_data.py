import pandas as pd
import argparse, re, os, vcf

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Process variation data from a file to produce csv of SNPs frequencies by country')
    parser.add_argument('--input_path', '-i', help='path to the variation data file in vcf format', required=True)
    parser.add_argument('--output_path', '-o',
                        help='path to the output collecting SNP frequencies by countries in csv format', required=True)
    parser.add_argument('--metadata_to_group_by', '-f',
                        help='metadata field to group by the SNPs frequencies by. Options: region, country, age, sex. default is country',
                        required=False, default='country')
    parser.add_argument('--metadata_path', '-m', help='path ro the metadata file', required=True)

    args = parser.parse_args()
    input_path = args.input_path
    output_path = args.output_path
    metadata_to_group_by = args.metadata_to_group_by
    metadata_path = args.metadata_path

    # read input file
    vcf_reader = vcf.Reader(open(input_path))
    records = [record for record in vcf_reader]

    metadata = pd.read_csv(metadata_path, sep='\t')
    if metadata_to_group_by == 'division':
        metadata = metadata.loc[(metadata['country'] == 'USA')]
        print("filtered df size: ", metadata.shape)
    samples_num = metadata.shape[0] # this is the number of samples, including those that demonstrated complete identity with the reference genome

    # get the count of each metadata value across samples
    metadata_to_samples_count = dict()
    metadata_values = metadata[metadata_to_group_by].unique()
    for metadata_value in metadata_values:
        metadata_to_samples_count[metadata_value] = metadata.loc[metadata[metadata_to_group_by] == metadata_value].shape[0]

    # collect the SNPs and their frequencies across all the samples
    SNPs_to_frequency = dict()
    for record in records:
        if (str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT)) not in SNPs_to_frequency:
            SNPs_to_frequency[(str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT))] = 1
        else:
            SNPs_to_frequency[(str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT))] += 1
    for SNP in SNPs_to_frequency:
        SNPs_to_frequency[SNP] /= samples_num


    # search for SNPs with at least 5% frequency across worldwide population
    frequent_SNPs = []
    for SNP in SNPs_to_frequency:
        if SNPs_to_frequency[SNP] > 0.01:
            frequent_SNPs.append(SNP)

    # now, collect a the frequencies of the filtered SNPs across countries into a dataframe
    data = dict()
    for record in records:
        sample_id_regex = re.compile("_sample_id_(.*?)_date_", re.DOTALL)
        sample_id = sample_id_regex.search(record.INFO['QNAME']).group(1)
        try:
            metadata_value = metadata.loc[(metadata['gisaid_epi_isl'] == sample_id)][metadata_to_group_by].item()
            if not metadata_value in data:
                data[metadata_value] = dict()
                for SNP in frequent_SNPs:
                    data[metadata_value][SNP] = 0
            if (str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT)) in frequent_SNPs:
                data[metadata_value][(str(record.CHROM), str(record.POS), str(record.REF), str(record.ALT))] += 1
        except:
            print("sample " + sample_id + " is not from USA, so it is being skipped")
            continue
    for metadata_value in metadata_to_samples_count:
        for SNP in frequent_SNPs:
            try:
                data[metadata_value][SNP] /= metadata_to_samples_count[metadata_value]
            except:
                print("no samples with variants for ", metadata_to_group_by, " : ", metadata_value)
                data[metadata_value] = dict()
                for SNP in frequent_SNPs:
                    data[metadata_value][SNP] = 0

    # insert data to a dataframe
    df = pd.DataFrame(columns=[metadata_to_group_by, 'frequency_across_samples'] + [
        'chrom_' + str(frequent_SNP[0]) + '_pos_' + frequent_SNP[1] + '_ref_' + frequent_SNP[2] + '_alt_' +
        frequent_SNP[3] for frequent_SNP in frequent_SNPs])
    for metadata_value in metadata_to_samples_count:
        value = {metadata_to_group_by: metadata_value}
        value['frequency_across_samples'] = metadata_to_samples_count[metadata_value] / samples_num
        for SNP in frequent_SNPs:
            value['chrom_' + SNP[0] + '_pos_' + SNP[1] + '_ref_' + SNP[2] + '_alt_' + SNP[3]] = data[metadata_value][
                SNP]
        df = df.append(value, ignore_index=True)
    df.to_csv(output_path)