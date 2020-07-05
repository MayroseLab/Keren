# Imports
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from lmfit import Model
from scipy.stats import pearsonr
import time

COUNTRY = "Hong Kong"
GEOTYPE = "sub-region" # country/region | city | sub-region
TRANSPORTATION_TYPE = "driving" # walking | driving | transit
METHOD = 0
MULTIPLESP = False
DEBUG = False
# OECD_countries = ['Austria', 'Belgium', 'Chile', 'Colombia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Israel', 'Italy', 'Japan', 'Latvia', 'Lithuania', 'Luxembourg', 'Mexico', 'Netherlands', 'New Zealand', 'Norway', 'Poland', 'Portugal', 'Slovenia', 'Slovakia', 'South Korea', 'Spain', 'Sweden', 'Switzerland', 'Turkey', 'United Kingdom', 'Australia', 'Canada', 'United States']
US_states = ['Alabama','Alaska','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Florida','Georgia','Hawaii','Idaho','Illinois', 'Indiana','Iowa','Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi','Missouri','Montana', 'Nebraska','Nevada','New Hampshire','New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio','Oregon','Pennsylvania', 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming', 'Oklahoma']


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
# the day_of_10_deaths is given relative to start date of death data collecting (22.1) while the mobility dat starts at 13.1
# the fearutes representing an absolute time (e.g., begin_quarantine) need to be standardized to the day_of_10_deaths
def extract_features_per_state(df, day_of_10_deaths):
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()
    L, x0, k, b, x1, a = fit_to_sigplus(df)
    try:
        f_x1 = sigmoid_plus(xdata,L, x0, k, b, x1, a)[int(x1)]
    except:
        f_x1 = b
    day_of_10_deaths_in_mobility_scale = day_of_10_deaths + 9                # Apple data starts at 13/1, while death data starts at 22/1, so 9 days difference
    #features = {"Day entering lockdown (shifted by first 10 deaths)": x0 + np.log(19)/k - day_of_10_deaths_in_mobility_scale,    # np.log(19)=ln(0.95/0.05), scale absolute date value by day_of_10_deaths_in_mobility_scale
    features = {"Day entering lockdown": x0 + np.log(19) / k,
                #"Day lockdown was fully invoked (shifted by first 10 deaths)": x0 - np.log(19)/k - day_of_10_deaths_in_mobility_scale,    # scale absolute date value by day_of_10_deaths_in_mobility_scale
                "Day lockdown was fully invoked": x0 - np.log(19) / k,
                "# days for lockdown to be invoked": abs(2*np.log(19) / k),
                "Lockdown magnitude": L/(L+f_x1), # originaly was L/(L+b), altered to assure sigmoid full description in the sigmoid+
                "# days under lockdown": x1 - (x0 - np.log(19)/k),
                "t1": x1, #- day_of_10_deaths_in_mobility_scale,        # scale absolute date value by day_of_10_deaths_in_mobility_scale
                "a": a}
    # features["increase in mobility before lockdown"] = (df[date_cols[int(features["Day entering lockdown (shifted by first 10 deaths)"]) - 14: int(features["Day entering lockdown (shifted by first 10 deaths)"]) - 7]].mean(axis=1) - 100).values[0]
    # features["variability in mobility during lockdown"] = np.var(df[date_cols[int(features["Day lockdown was fully invoked (shifted by first 10 deaths)"]):int(features["t1"])]].values)
    features["increase in mobility before lockdown"] = (df[date_cols[int(features["Day entering lockdown"]) - 14: int(features["Day entering lockdown"]) - 7]].mean(axis=1) - 100).values[0]
    features["variability in mobility during lockdown"] = np.var(df[date_cols[int(features["Day lockdown was fully invoked"]):int(features["t1"])]].values)
    pearson, pvalue = pearsonr(ydata,sigmoid_plus(xdata,L, x0, k, b, x1, a))
    features["pearson_of_fit *^ 2"] = pearson**2
    features["pvalue_of_fit"] = pvalue
    if np.isnan(features["variability in mobility during lockdown"]):
        features["variability in mobility during lockdown"] = 0
    return features

def fit_to_sigplus(df, method = METHOD, multiple_starting_points=MULTIPLESP):

    # extract data to fit the function to
    date_cols = [x for x in df.columns if "2020" in x]
    xdata = np.arange(len(date_cols))
    ydata = df[date_cols].values.flatten()

    # set starting points (options are 1 or multiple - 3)
    L_0 = (max(ydata) - min(ydata))
    x0_0 = np.argwhere(ydata < (100 + min(ydata))/2)[0][0]
    k_0 = -1
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
    L, x0, k, b, x1, a = fit_to_sigplus(filtered_df, method=method, multiple_starting_points=multiple_starting_points)
    cost = get_solution_cost(xdata, ydata, [L, x0, k, b, x1, a], cost_method=cost_method)
    return cost


if __name__ == "__main__":

    df = pd.read_csv("mobility_data/apple_data/applemobilitytrends-2020-05-10.csv")
    # change countries names to fit the rest of data sources: UK->United Kingdom, Republic of Korea->South Korea
    df['region'] = df['region'].replace({'UK': 'United Kingdom', 'Republic of Korea': 'South Korea'})
    df = df.loc[df.geo_type == GEOTYPE]
    # death_data = pd.read_excel("../quarantine-minipaper/data csv files/table_match_to_death_sigmoid_v2.xlsx")


    if DEBUG:
        columns = ['country', 'cost_of_1sp_old', 'cost_of_3sp_old', 'cost_of_1sp_new', 'cost_of_3sp_new', 'time_to_eval']
        categories = df.region.unique()
        costs_dict = dict()
        for category in categories:
            if type(category) == str:
                print("category\t", category)
                start = time.time()
                filtered_df = df[(df.region == category) & (df.transportation_type == TRANSPORTATION_TYPE)]
                try:
                    cost_of_1sp_old = evaluate_fit(filtered_df, method=0, multiple_starting_points=False)
                except Exception as e:
                    print('category: ', category, ', error in using curve_fit with one starting point: ', e)
                    cost_of_1sp_old = float('nan')
                try:
                    cost_of_3sp_old = evaluate_fit(filtered_df, method=0, multiple_starting_points=True)
                except Exception as e:
                    print('category: ', category, ', error in using curve_fit with 3 starting points: ', e)
                    cost_of_3sp_old = float('nan')
                try:
                    cost_of_1sp_new = evaluate_fit(filtered_df, method=1, multiple_starting_points=False)
                except Exception as e:
                    print('category: ', category, ', error in using spmodel.fit with one starting point: ', e)
                    cost_of_1sp_new = float('nan')
                try:
                    cost_of_3sp_new = evaluate_fit(filtered_df, method=1, multiple_starting_points=True)
                except Exception as e:
                    print('category: ', category, ', error in using spmodel.fit with 3 starting points: ', e)
                    cost_of_3sp_new = float('nan')
                end = time.time()
                costs_dict[category] = {'cost_of_1sp_old': cost_of_1sp_old, 'cost_of_3sp_old': cost_of_3sp_old, 'cost_of_1sp_new': cost_of_1sp_new, 'cost_of_3sp_new': cost_of_3sp_new, 'time_to_eval': end-start}
        cost_by_country = pd.DataFrame(costs_dict).T
        cost_by_country.to_csv("mobility_data/fit_eval_{}-2020-05-18.csv".format(GEOTYPE.replace("/","_or_") + "_" + TRANSPORTATION_TYPE))

    else:
        columns = [ "country",
                    # "Day entering lockdown (shifted by first 10 deaths)",                     # float representing the absolute start time of quarantine
                    "Day entering lockdown",
                    # "Day lockdown was fully invoked (shifted by first 10 deaths)",                      # float representing the absolute time in which the quarantine has reached maximal restriction of movement
                    "Day lockdown was fully invoked"
                    "# days for lockdown to be invoked",              # float representing the absolute time it took from when the qurantine was applied to maximal enforcement
                    "Lockdown magnitude",               # float representing the magnitude of drop in the mobility during the enforcement of the quarantine
                    "# days under lockdown",            # float representing the duration in which the minimal mobility obtained from full quarantine was maintained
                    "t1",                  # float representing the absolute time in which the mobility started to increase after reaching its minimal value
                    "a",   # float representing the rate of increase in mobility following the end of the quarantine
                    "increase in mobility before lockdown",                 # float representing the increase in mobility before the quarantine has begun
                    "variability in mobility during lockdown",               # variation of the mobility during the time in which the minimal mobility obtained from full quarantine was maintained
                    "pearson_of_fit^2",                       # measure of the quality of fit of the sigmoidplus function to the data. Should be used to measure the credibility of the feature values of the respective categories
                    "pvalue_of_fit"]                        # measure of the significance of the quality of fit of the sigmoidplus function to the data. Should be used to measure the credibility of the feature values of the respective categories
        features_dict = {}
        categories = US_states # OECD_countries
        # print("category\tfailure_cause")
        for category in categories:
            if type(category) == str:
                print("category: ", category)
                # try:
                filtered_df = df[(df.region == category) & (df.transportation_type == TRANSPORTATION_TYPE)]
                day_of_10_deaths = float("nan") # death_data.loc[death_data['Country name and its population'] == category]['Calculated directly from death data w.o sigmoid model'].values[0]
                features_dict[category] = extract_features_per_state(filtered_df, day_of_10_deaths)
                # except Exception as e:
                #     print(category, "\t", e)

        features_by_country = pd.DataFrame(features_dict).T
        features_by_country.index.names = ['country']
        features_by_country.to_csv("mobility_data/us_states_features_{}-2020-05-10.csv".format(
            GEOTYPE.replace("/", "_or_") + "_" + TRANSPORTATION_TYPE))
