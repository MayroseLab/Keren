# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy.signal import medfilt


COUNTRY = "Belgium"
GEOTYPE = "city" # country/region | city | sub-region
TRANSPORTATION_TYPE = "driving" # walking | driving | transit

def sigmoid(x, L, x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y


def sigmoid_plus(x, L, x0, k, b, x1, a):
    """
    :param x:
    :param L: start of sigmoid (on y-axis)
    :param x0: the middle point of decline
    :param k: 2/k is the range of decline
    :param x1: end of sigmoid (point on x-axis)
    :param a: linear function: the slope of incline following sigmoid drop
    :param b: sigmoid+linear function: function shift along y axis
    :return:
    """
    y1 = L / (1 + np.exp(-k * (x - x0))) + b
    y2 = a * x - a * x1
    y = y1
    if (int(x1) + 1) < len(x) and (x1 > (x0 - 1 / k)):
        y[(int(x1) + 1):] = y[(int(x1) + 1):] + y2[(int(x1) + 1):]
    return y


def extract_features_per_state(filtered_df):
    date_cols = [x for x in filtered_df.columns if "2020" in x]
    L, x0, k, b, x1, a = fit_to_sigplus(filtered_df)
    features = {"begin_quarantine": x0 + 1 / k,
                 "full_quarantine": x0 - 1 / k,
                 "time_to_full_quarantine": -2 / k,
                 "strength_of_quarantine": L,
                 "length_of_full_quarantine": x1 - (x0 - 1 / k),
                 "day_quit_quarantine": x1,
                 "increase_rate_from_quit_quarantine": a}
    features["increase_before_drop"] = (filtered_df[date_cols[int(features["begin_quarantine"]) - 14: int(features["begin_quarantine"]) - 7]].mean(axis=1) - 100).values[0]
    features["quarantine_enforcement"] = np.var(filtered_df[date_cols[int(features["full_quarantine"]):int(features["day_quit_quarantine"])]].values)
    if np.isnan(features["quarantine_enforcement"]):
        features["quarantine_enforcement"] = 0

    return features


def fit_to_sigplus(df):
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()
    try:
        L_0 = (max(ydata) + min(ydata)) / 2
    except Exception as e:
        print("e: ", e)
        print("ydata: ", ydata)
        print("date_cols: ", date_cols)
    x0_0 = np.argwhere(ydata < (100 + min(ydata)) / 2)[0][0]
    x1_0 = len(date_cols) - 30  # x0_0 + np.argwhere(ydata_tmp < (100+min(ydata_tmp))/2)[0][0]
    a_0 = 1
    k_0 = -2
    b_0 = (min(ydata))
    popt, _ = curve_fit(sigmoid_plus, xdata, ydata, [L_0, x0_0, k_0, b_0, x1_0, a_0], method="dogbox", maxfev=10000)
    (L, x0, k, b, x1, a) = popt
    # # I THINK THIS SHOULD BE PLACED BACK IN SIGMOID_PLUS AND APPLY THE ASSERTION ON ANY CHOSEN SET OF PARAMETERS
    # # verify credibility of optimized function parameters
    # assert (x0 > 0), "inferred beginning of quarantine does not exist within date range"
    # assert (L > 0), "inferred strength of quarantine indicates increase in mobility rather than decrease"
    # assert (k < 0), "inferred direction of sigmoid indicates increase in mobility rather than decrease"
    # assert (x1 > (x0 - 1 / k)), "inferred parameters suggest the the quarantine has ended before it has begun"

    return L, x0, k, b, x1, a


if __name__ == "__main__":

    df = pd.read_csv("mobility_data/applemobilitytrends-2020-05-18.csv")
    df.drop(columns=["11/05/2020", "12/05/2020"], inplace=True) # drop dates with missing data across all countries, as indicated by https://www.apple.com/covid19/mobility
    df = df.loc[df.geo_type == GEOTYPE]
    columns = [ "country",
                "begin_quarantine",                     # float representing the absolute start time of quarantine
                "full_quarantine",                      # float representing the absolute time in which the quarantine has reached maximal restriction of movement
                "time_to_full_quarantine",              # float representing the absolute time it took from when the qurantine was applied to maximal enforcement
                "strength_of_quarantine",               # float representing the magnitude of drop in the mobility during the enforcement of the quarantine
                "length_of_full_quarantine",            # float representing the duration in which the minimal mobility obtained from full quarantine was maintained
                "day_quit_quarantine",                  # float representing the absolute time in which the mobility started to increase after reaching its minimal value
                "increase_rate_from_quit_quarantine",   # float representing the rate of increase in mobility following the end of the quarantine
                "increase_before_drop",                 # float representing the increase in mobility before the quarantine has begun
                "quarantine_enforcement"]               # variation of the mobility during the time in which the minimal mobility obtained from full quarantine was maintained
    features_dict = {}
    categories = df.region.unique()

    for category in categories:
        if category in ['Belo Horizonte', 'Aachen']: # Singapore - failed to optimize parameters,
            continue
        if type(category) == str:
            print("category: ", category)
            filtered_df = df[(df.region == category) & (df.transportation_type == TRANSPORTATION_TYPE)]
            features_dict[category] = extract_features_per_state(filtered_df)

    features_by_country = pd.DataFrame(features_dict).T
    features_by_country.to_csv("mobility_data/{}-2020-05-18.csv".format(GEOTYPE.replace("/","_or_") + "_" + TRANSPORTATION_TYPE))