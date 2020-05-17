import pandas as pd
import argparse, re, os

if __name__ == '__main__':

    os.chdir(r"C:\Users\ItayMNB7\Google Drive\PhD\ML_Covid19\GISAID")

    # process input from command line
    parser = argparse.ArgumentParser(description='Process variation data from a file to produce csv of SNPs frequencies by country')
    parser.add_argument('--input_path', '-i', help='path to the variation data file', required=False, default="variation_data.xlsx")
    parser.add_argument('--variants_to_collect', '-a', help='path to a vcf file with the variants that data should be collected for', required=False, default="frequent_variants.tsv")
    parser.add_argument('--output_path', '-o', help='path to the CSV output', required=False, default="variants_frequencies_by_country.csv")

    args = parser.parse_args()
    input_path = args.input_path
    variants_to_collect = args.variants_to_collect
    output_path = args.output_path

    # read input file
    df = pd.read_excel(input_path)

    # read varaints of interest
    variants = []
    with open(variants_to_collect, "r") as infile:
        lines = infile.readlines()
    for line in lines:
        components = line.split("\t")
        variants.append((components[0], components[1], components[2].replace("\n", "")))
    print("variants: ", variants)

    # create a new dataframe, in which each entry depicts the frequencies of the 20 variants (named by their positions) across the country
    countries = []
    country_regex = re.compile("hCoV-19_(.*?)[_|;]", re.IGNORECASE | re.MULTILINE | re.DOTALL)
    for index, row in df.iterrows():
        info = str(row[['INFO']])
        country = country_regex.search(info).group(1)
        countries.append(country)
    df.insert(2, "Country", countries, True)
    unique_countries = df.Country.unique()

    # count the number of records per country and position
    columns = ['Country'] + [variant[0] + ":" + variant[1] + "->" + variant[2] for variant in variants]
    variants_frequencies = pd.DataFrame(columns=columns)
    for country in unique_countries:
        record = {'Country': country}
        for variant in variants:
            filtered_df = df.loc[(df['Country'] == country) & (df['POS'] == variant[0])]
            print(filtered_df.shape[0])
            frequency = (df.loc[(df['Country'] == country) & (df['POS'] == variant[0]) & (df['REF'] == variant[1]) & (df['ALT'] == variant[2])]).shape[0]
            record[variant[0] + ":" + variant[1] + "->" + variant[2]] = frequency
        variants_frequencies = variants_frequencies.append(record, ignore_index=True)
    variants_frequencies.to_csv(output_path)







