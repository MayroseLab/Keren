# Imports
import pandas as pd
from scipy.stats import pearsonr

GEOTYPE = "country/region" # country/region | city | sub-region
TRANSPORTATION_TYPE_1 = "walking" # walking | driving | transit
TRANSPORTATION_TYPE_2 = "driving" # walking | driving | transit
OECD_countries = ['Austria', 'Belgium', 'Chile', 'Colombia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Israel', 'Italy', 'Japan', 'Latvia', 'Lithuania', 'Luxembourg', 'Mexico', 'Netherlands', 'New Zealand', 'Norway', 'Poland', 'Portugal', 'Slovenia', 'Slovakia', 'South Korea', 'Spain', 'Sweden', 'Switzerland', 'Turkey', 'United Kingdom', 'Australia', 'Canada', 'United States']

if __name__ == "__main__":

    df = pd.read_csv("mobility_data/applemobilitytrends-2020-05-10.csv")
    df['region'] = df['region'].replace({'UK': 'United Kingdom', 'Republic of Korea': 'South Korea'})
    df_1 = df.loc[(df['region'].isin(OECD_countries)) & (df.geo_type == GEOTYPE) & (df.transportation_type == TRANSPORTATION_TYPE_1)]
    df_2 = df.loc[(df['region'].isin(OECD_countries)) & (df.geo_type == GEOTYPE) & (df.transportation_type == TRANSPORTATION_TYPE_2)]

    subjects_for_correlation = df.region.unique() # [x for x in df.columns if "2020" in x]
    elements = dict()
    for subject in subjects_for_correlation:
        element = dict()
        pearson, pvalue = pearsonr(df_1[subject].values.flatten(), df_2[subject].values.flatten())
        element["pearson_correlation ^ 2 "] = pearson ** 2
        element["pvalue"] = pvalue
        elements[subject] = element


    corr_df = pd.DataFrame(elements).T
    # corr_df.index.names = [index]
    corr_df.to_csv("../quarantine-minipaper/correlations_of_mobility_data_{}_oecd.csv".format(TRANSPORTATION_TYPE_1 + "_" + TRANSPORTATION_TYPE_2))
