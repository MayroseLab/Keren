# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import inspect
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy.signal import medfilt
import datetime

GEOTYPE = 'country_region'  # country_region_code | country_region (country) | sub_region_1 (state) | sub_region_2 (county)
VALUE = 'workplaces_percent_change_from_baseline'   # retail_and_recreation_percent_change_from_baseline,
                                                    # | grocery_and_pharmacy_percent_change_from_baseline
                                                    # | parks_percent_change_from_baseline
                                                    # | transit_stations_percent_change_from_baseline
                                                    # | workplaces_percent_change_from_baseline
                                                    # | residential_percent_change_from_baseline

def reformat_data(input_path, output_dir, value=VALUE, group_by_value=GEOTYPE):
    df = pd.read_csv(input_path, header=0,
                     dtype={'country_region_code': str, 'country_region': str, 'sub_region_1': str, 'sub_region_2': str,
                            'date': str})
    # df = df[[group_by_value, 'date', value]].groupby(group_by_value, as_index=False).aggregate(np.sum)
    df = df[[group_by_value, 'date', value]].groupby([group_by_value, 'date'], as_index=False).aggregate(np.sum)
    out_df = df.pivot(index=group_by_value, columns='date', values=value)
    out_df.to_csv(output_dir + value + "_by_" + group_by_value + ".csv")
    return out_df

#### functions to fit data to ####

def sigmoid(x, L, x0, k, b):
    """
    :param x: x data
    :param L: start of sigmoid (on y-axis)
    :param x0: the middle point of decline
    :param k: 2/k is the range of decline
    :param b: sigmoid+linear function: function shift along y axis
    :return:
    """
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return y


def sigmoid_plus(x, L, x0, k, b, x1, a):
    """
    :param x: x data
    :param L: start of sigmoid (on y-axis)
    :param x0: the middle point of decline
    :param k: 2/k is the range of decline
    :param x1: end of sigmoid (point on x-axis)
    :param a: linear function: the slope of incline following sigmoid drop
    :param b: sigmoid+linear function: function shift along y axis
    :return:
    """
    y1 = sigmoid(x, L, x0, k, b)
    y2 = a * x - a * x1
    y = y1 # start with just sigmoid part
    if (int(x1) + 1) < len(x) and (x1 > (x0 - 1 / k)):
        y[(int(x1) + 1):] = y[int(x1)] + y2[(int(x1) + 1):] # plug in the plus part
    return y


def minus_sigmoid(x, L, x0, k, b, x1, a):
    """
     :param x: x data
     :param a: linear function: the slope of decline following sigmoid drop
     :param x0: end of linear function (point on x-axis)
     :param L: start of sigmoid (on y-axis)
     :param x1: the middle point of incline
     :param k: 2/k is the range of incline
     :param b: sigmoid+linear function: function shift along y axis
     :return:
     """
    y1 = -a * x[:int(x0)+1] - (-a) * x0
    y2 = sigmoid(x[(int(x0) + 1):], L, x1, k, b)
    y = y1  # start with just minus part
    if (int(x0) + 1) < len(x) and (x0 > (x1 - 1 / k)):
        y[(int(x0) + 1):] = y[int(x0)+1] + y2  # plug in the plus part
    return y


def plus_sigmoid(x, L, x0, k, b, x1, a):
    """
    :param x: x data
    :param a: linear function: the slope of incline following sigmoid drop
    :param x0: end of linear function (point on x-axis)
    :param L: start of sigmoid (on y-axis)
    :param x1: the middle point of incline
    :param k: 2/k is the range of incline
    :param b: sigmoid+linear function: function shift along y axis
    :return:
    """
    y1 = a * x[:int(x0)+1] - a * x0
    y2 = sigmoid(x[(int(x0) + 1):], L, x1, k, b)
    y = y1  # start with just minus part
    if (int(x0) + 1) < len(x) and (x0 > (x1 - 1 / k)):
        y[(int(x0) + 1):] = y[int(x0)] + y2  # plug in the plus part
    return y


value_to_function = {'retail_and_recreation_percent_change_from_baseline': minus_sigmoid,
                     'grocery_and_pharmacy_percent_change_from_baseline': plus_sigmoid,
                      'parks_percent_change_from_baseline': minus_sigmoid,
                      'transit_stations_percent_change_from_baseline': minus_sigmoid,
                      'workplaces_percent_change_from_baseline': minus_sigmoid,
                      'residential_percent_change_from_baseline': plus_sigmoid}

#### fit function to data ####

def set_starting_point(xdata, ydata, value):
    function = value_to_function[value]
    function_arguments = inspect.getfullargspec(function)
    L_0 = (max(ydata) + min(ydata)) / 2
    k_0 = 2
    b_0 = (min(ydata))

    starting_point = [xdata, L_0, 0, k_0, b_0, 0, 0]
    # fit to minus_sigmoid
    if value in ['retail_and_recreation_percent_change_from_baseline', 'parks_percent_change_from_baseline', 'transit_stations_percent_change_from_baseline', 'workplaces_percent_change_from_baseline']:
        x0_0 = np.argwhere(ydata < (100 + min(ydata)) / 2)[0][0]
        x1_0 = len(xdata) - 30  # x0_0 + np.argwhere(ydata_tmp < (100+min(ydata_tmp))/2)[0][0]
        a_0 = 1
        starting_point = [xdata, L_0, x0_0, k_0, b_0, x1_0, a_0]

    # fit to plus_sigmoid
    elif value in ['grocery_and_pharmacy_percent_change_from_baseline', 'residential_percent_change_from_baseline']:
        L_0 = (max(ydata) + min(ydata)) / 2
        x0_0 = np.argwhere(ydata < (100 + min(ydata)) / 2)[0][0]
        x1_0 = len(xdata) - 30  # x0_0 + np.argwhere(ydata_tmp < (100+min(ydata_tmp))/2)[0][0]
        a_0 = 1
        k_0 = -2
        b_0 = (min(ydata))
        starting_point = [xdata, L_0, x0_0, k_0, b_0, x1_0, a_0]
    return starting_point



def fit_to_function(df, value=VALUE): # the value should dictate the choice of starting point
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()
    function = value_to_function[value]
    strating_point = set_starting_point(xdata, ydata, value)
    popt, _ = curve_fit(function, xdata, ydata, strating_point, method="dogbox", maxfev=10000)
    (L, x0, k, b, x1, a) = popt
    # # I THINK THIS SHOULD BE PLACED BACK IN SIGMOID_PLUS AND APPLY THE ASSERTION ON ANY CHOSEN SET OF PARAMETERS
    # # verify credibility of optimized function parameters
    # assert (x0 > 0), "inferred beginning of quarantine does not exist within date range"
    # assert (L > 0), "inferred strength of quarantine indicates increase in mobility rather than decrease"
    # assert (k < 0), "inferred direction of sigmoid indicates increase in mobility rather than decrease"
    # assert (x1 > (x0 - 1 / k)), "inferred parameters suggest the the quarantine has ended before it has begun"

    return L, x0, k, b, x1, a

if __name__ == '__main__':

    input_path = "mobility_data/google_Global_Mobility_Report_19-5-20.csv"
    output_dir = "mobility_data/"
    df = reformat_data(input_path, output_dir, value=VALUE, group_by_value=GEOTYPE)
