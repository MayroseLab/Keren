import re, pandas as pd, os, sys, argparse
from scipy.stats.distributions import chi2

def doLRT(null_logl, alternative_logl, df=1):
    LR = 2 * (alternative_logl - null_logl)
    pvalue = chi2.sf(LR, df)
    return pvalue


# integrate into main
def extract_hyphy_data(input_dir):

    colname_to_regex = {"dataset_id": re.compile("(COG.*?)\.fa", re.MULTILINE | re.DOTALL),
                        "null_logl": re.compile("Fitting the null \(K \:= 1\) model.*?Log\(L\)\s*=\s*(-\d*\.?\d*)",
                                                re.MULTILINE | re.DOTALL),
                        "alternative_logl": re.compile("Fitting the alternative model.*?Log\(L\)\s*=\s*(-\d*\.?\d*)",
                                                       re.MULTILINE | re.DOTALL),
                        "pvalue": re.compile("Likelihood ratio test\s*\*\*p\s*=\s*(\d*\.?\d*)\*\*",
                                             re.MULTILINE | re.DOTALL),
                        "null_omega0": re.compile("Fitting the null(.*?\|){12}\s*(\d*\.?\d*)",
                                                  re.MULTILINE | re.DOTALL),
                        "null_omega1": re.compile("Fitting the null(.*?\|){17}\s*(\d*\.?\d*)",
                                                  re.MULTILINE | re.DOTALL),
                        "null_omega2": re.compile("Fitting the null(.*?\|){22}\s*(\d*\.?\d*)",
                                                  re.MULTILINE | re.DOTALL),
                        "null_p0": re.compile("Fitting the null(.*?\|){13}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                        "null_p1": re.compile("Fitting the null(.*?\|){18}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                        "alternative_omega0": re.compile(
                            "Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){12}\s*(\d*\.?\d*)",
                            re.MULTILINE | re.DOTALL),
                        "alternative_omega1": re.compile(
                            "Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){17}\s*(\d*\.?\d*)",
                            re.MULTILINE | re.DOTALL),
                        "alternative_omega2": re.compile(
                            "Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){22}\s*(\d*\.?\d*)",
                            re.MULTILINE | re.DOTALL),
                        "alternative_p0": re.compile(
                            "Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){13}\s*(\d*\.?\d*)",
                            re.MULTILINE | re.DOTALL),
                        "alternative_p1": re.compile(
                            "Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){18}\s*(\d*\.?\d*)",
                            re.MULTILINE | re.DOTALL),
                        "alternative_k": re.compile("Relaxation\/intensification parameter\s*\(K\)\s*=\s*(\d*\.?\d*)",
                                                    re.MULTILINE | re.DOTALL),
                        "job_id": re.compile("(\d*)\.power8", re.MULTILINE | re.DOTALL)}
    colnames = list(colname_to_regex.keys()) + ["LR", "significant", "job_id"]
    df = pd.DataFrame(colnames)

    for path in os.listdir(input_dir):
        record = dict()
        record["job_id"] = colname_to_regex["job_id"].search(path).group(1)
        with open(input_dir+path, "r") as infile:
            content = infile.read()
        for colname in list(colname_to_regex.keys()):
            # print("colname: ", colname, "\nregex: ", colname_to_regex[colname].pattern, "\npath: ", input_dir+path)
            if "omega" in colname:
                record[colname] = float(colname_to_regex[colname].search(content).group(2))
            elif "p0" in colname or "p1" in colname:
                record[colname] = float(colname_to_regex[colname].search(content).group(2)) / 100
            elif colname != "job_id" and colname != "dataset_id":
                record[colname] = float(colname_to_regex[colname].search(content).group(1))
            elif colname == "dataset_id":
                record[colname] = colname_to_regex[colname].search(content).group(1)
        record["LR"] = 2 * (record["alternative_logl"] - record["null_logl"])
        pvalue = doLRT(record["null_logl"], record["alternative_logl"])
        if abs(pvalue - record["pvalue"]) > 0.05:
            print("error! coputed pvalue is different than the one reported by hyhphy!")
            print("input path: ", input_dir + path)
            print("replicate: ", record["replicate"])
            print("computed pvalue: ", pvalue)
            print("reported pvalue: ", record["pvalue"])
            exit(1)
        if pvalue < 0.05:
            record["significant"] = 1
        else:
            record["significant"] = 0
        df = df.append(record, ignore_index=True)
    return df

if __name__ == '__main__':
    # process input from command line
    parser = argparse.ArgumentParser(
        description='Extract the results of RELAX analysis in Bio++ on real datasets')
    parser.add_argument('--input_dir', '-i', help='directory that holds the output files by relax to be analyzed',
                        required=True)
    parser.add_argument('--output_dir', '-o', help='directory that will hold the csv output of the analysis',
                        required=True)

    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir

    # process data
    df = extract_hyphy_data(input_dir)

    df.to_csv(output_dir + "res.csv")

