import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# pd.read_csv(measurements_path)


# plot vector x against y according to a vector x and a vector y of values, and a list of relevant csv paths, and the function required (mean/median ect.)
def plot_data(x_vector, input_files_prefix, input_files_suffix, plot_title, x_title, y_title, y_name, output_path, func="mean"):
    y_vector = np.array()
    sd_vector = np.array()
    for x in x_vector:
        input_path = input_files_prefix + str(x) + input_files_suffix
        df = pd.read_csv(input_path, sep="\t")
        values = np.array(df[y_name])
        func_val = 0
        if func == "mean":
            func_val = np.mean(values)
        elif func == "median":
            func_val = np.median(values)
        sd_vector.append(np.std(values))
        y_vector.append(func_val)
    plt.plot(x_vector, y_vector)
    plt.title(plot_title)
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    for sd in sd_vector: # add standard deviation lines to the plot for each y_val
        plt.axvline(x=sd, color='k', linestyle='--')
    plt.savefig(output_path)
    return 0








