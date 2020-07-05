# Imports
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from lmfit import Model
from scipy.stats import pearsonr
import time

OECD_members = ['Austria', 'Belgium', 'Canada', 'Chile', 'Colombia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Israel', 'Italy', 'Japan', 'Republic of Korea', 'Latvia', 'Lithuania', 'Luxembourg', 'Mexico', 'Netherlands', 'New Zealand', 'Norway', 'Poland', 'Portugal', 'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Turkey', 'UK', 'United States']
GEOTYPE = "country" # country/region | city | sub-region
TRANSPORTATION_TYPE = "driving" # walking | driving | transit
METHOD = 0
MULTIPLESP = False


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
    sig_num = 3
    y1 = sigmoid(x, L, x0, k, b)
    if (int(x1) + 1) < len(x) and (x1>(x0-sig_num/k)):
        y2 = a * x - a * x1 + y1[int(x1)] # y1[int(x1)] is the y point from which the sigmoid ends
    y = y1
    if (int(x1) + 1) < len(x) and (x1>(x0-sig_num/k)):
        y[(int(x1) + 1):] = y2[(int(x1) + 1):]
    return y

# extracts features based on sigmoid plus fit
def extract_features_per_state(df):
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()
    L, x0, k, b, x1, a = fit_to_sigplus(df)
    features = {"begin_quarantine": x0 + 3 / k,
                 "full_quarantine": x0 - 3 / k,
                 "time_to_full_quarantine": abs(6 / k),
                 "strength_of_quarantine": L/(L+b),
                 "length_of_full_quarantine": x1 - (x0 - 3 / k),
                 "day_quit_quarantine": x1,
                 "increase_rate_from_quit_quarantine": a}
    features["increase_before_drop"] = (df[date_cols[int(features["begin_quarantine"]) - 14: int(features["begin_quarantine"]) - 7]].mean(axis=1) - 100).values[0]
    features["quarantine_enforcement"] = np.var(df[date_cols[int(features["full_quarantine"]):int(features["day_quit_quarantine"])]].values)
    features["pearson_of_fit"] = pearsonr(ydata,sigmoid_plus(xdata,L, x0, k, b, x1, a))[0]
    if np.isnan(features["quarantine_enforcement"]):
        features["quarantine_enforcement"] = 0
    return features

def fit_to_sigplus(df, method = METHOD, multiple_starting_points=MULTIPLESP):

    # extract data to fit the function to
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()

    # set starting points (options are 1 or multiple - 3)
    L_0 = (max(ydata) - min(ydata))
    x0_0 = np.argwhere(ydata < (100 + min(ydata))/2)[0][0]
    k_0 =  -3 # because sigmoid(x,k) = 0.95 when kx = -3
    b_0 = (min(ydata))
    x1_0 = len(date_cols) - 30 # assumes sigmoid ends at the end of the data (i.e., quarantine wasn't over as soon as it reached maximal capacity)
    a_0 = 1
    starting_points = [[L_0, x0_0, k_0, b_0, x1_0, a_0]]

    if multiple_starting_points:
        L_1 = max(ydata)
        x0_1 = np.argwhere(ydata == np.median(ydata))[0][0]
        k_1 = k_0
        b_1 = b_0
        x1_1 = x1_0
        a_1 = 0.5
        starting_points.append([L_1, x0_1, k_1, b_1, x1_1, a_1])
        L_2 = min(ydata)
        x0_2 = x0_1
        k_2 = k_0
        b_2 = b_0
        x1_2 = x1_0
        a_2 = 2
        starting_points.append([L_2, x0_2, k_2, b_2, x1_2, a_2])

    # fit the function
    solutions = dict()
    for i in range(len(starting_points)):
        if method == 0: # uses non linear least squares via function curve_fit with lm method
            try:
                solution, _ = curve_fit(sigmoid_plus, xdata, ydata, starting_points[i], method="lm", maxfev=100000)
                solution = tuple(solution)
            except:
                continue
        else:
            try:
                spmodel = Model(sigmoid_plus)
                params = spmodel.make_params(L=starting_points[i][0], x0=starting_points[i][1], k=starting_points[i][2], b=starting_points[i][3], x1=starting_points[i][4], a=starting_points[i][5])
                result = spmodel.fit(ydata, params, x=xdata)
                solution = (result.best_values['L'], result.best_values['x0'], result.best_values['k'], result.best_values['b'], \
                       result.best_values['x1'], result.best_values['a'])
            except:
                continue
        cost = get_solution_cost(xdata, ydata, solution)
        solutions[solution] = cost
    optimal_solution = min(solutions, key=solutions.get) # take the optimal solution to be the one with the minimal cost
    L_opt = optimal_solution[0]
    x0_opt = optimal_solution[1]
    k_opt = optimal_solution[2]
    b_opt = optimal_solution[3]
    x1_opt = optimal_solution[4]
    a_opt = optimal_solution[5]
    # plt.figure()
    # plt.plot(xdata, ydata, 'b')
    # plt.plot(xdata, sigmoid_plus(xdata, L_opt, x0_opt, k_opt, b_opt, x1_opt, a_opt), 'r')
    # plt.plot([x0_opt + 2 / k_opt] * 200, range(0, 200), 'g')
    # plt.plot([x0_opt - 2 / k_opt] * 200, range(0, 200), 'y')
    # plt.plot([x0_opt] * 200, range(0, 200), 'p')
    # plt.show()
    return L_opt, x0_opt, k_opt, b_opt, x1_opt, a_opt

# copy from least_squares.py
def _wrap_func(func, xdata, ydata):
    def func_wrapped(params):
        return func(xdata, *params) - ydata
    return func_wrapped

# compute the least squares based cost of a solution
# cost methods: 0 - least square, 1 - 1-pearson
def get_solution_cost(x, y, solution, func=sigmoid_plus, cost_method=0):
    L = solution[0]
    x0 = solution[1]
    k = solution[2]
    b = solution[3]
    x1 = solution[4]
    a = solution[5]
    if cost_method == 0:
        fun = _wrap_func(func, x, y)
        np.atleast_1d(fun(solution))  # error here : operands could not be broadcast together with shapes (125,) (19125,)
        cost = 0.5 * np.dot(sigmoid_plus(x, L, x0, k, b, x1, a), sigmoid_plus(x, L, x0, k, b, x1, a))
    else:
        cost = 1-pearsonr(y, sigmoid_plus(x, L, x0, k, b, x1, a))[0]
    return cost

# evaluates the "cost" of the function in the optimal point
def evaluate_fit(df, method, multiple_starting_points, cost_method=1):
    value = dict()
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()
    L, x0, k, b, x1, a = fit_to_sigplus(df, method=method, multiple_starting_points=multiple_starting_points)
    cost = get_solution_cost(xdata, ydata, [L, x0, k, b, x1, a], cost_method=cost_method)
    return cost

if __name__ == '__main__':

    df = pd.read_csv("C:/Users/ItayMNB7/Dropbox/covid19/data/mobility_data/apple_features/features_country_or_region_{}-2020-05-10.csv".format(TRANSPORTATION_TYPE), dtype={'country': str})
    data = pd.read_csv("C:/Users/ItayMNB7/Dropbox/covid19/data/mobility_data/applemobilitytrends-2020-05-10.csv")
    df.sort_values(by=['pearson_of_fit'], inplace=True)  # sort from smallest (worst) pearson to best
    OECD_df = df.loc[df.country.isin(OECD_members)] # select OECD countries

    countries = ['Italy', 'Austria', 'Belgium', 'Denmark', 'France', 'Germany', 'Norway', 'Spain', 'Sweden', 'Switzerland', 'UK']
    fig, axis = plt.subplots(nrows=3, ncols=4, sharey='none', sharex='none', figsize=[5 * 12 + 2, 6 * 7.58 + 2], frameon=True)
    for i in range(11):
        ax = axis[int(i/4)][i % 4]
        #country = list(OECD_df.country)[i]
        country = countries[i]
        ax.set_title(country)
        country_data = data[(data.region == country) & (data.transportation_type == TRANSPORTATION_TYPE)]
        L_opt, x0_opt, k_opt, b_opt, x1_opt, a_opt = fit_to_sigplus(country_data)
        date_cols = [x for x in country_data.columns if "2020" in x]
        xdata = np.arange(len(date_cols))
        ydata = country_data[date_cols].values.flatten()
        joint_data = (country_data[date_cols]).T
        joint_data.columns = ['mobility']
        joint_data.index.names = ['date']
        joint_data = joint_data.reset_index()
        joint_data.sort_values(by=['mobility'], inplace=True, ascending=False)
        print('country: ', country, 'data:\n', joint_data.head())
        ax.plot(xdata, ydata, 'b')
        ax.plot(xdata, sigmoid_plus(xdata, L_opt, x0_opt, k_opt, b_opt, x1_opt, a_opt), 'r')
    plt.show()



