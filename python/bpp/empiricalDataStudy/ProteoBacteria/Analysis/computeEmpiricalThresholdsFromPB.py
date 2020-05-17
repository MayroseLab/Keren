import os, argparse, pandas as pd, numpy as np

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='extracts the empirical LR thresholds obtained by parametric bootstrapping and update them + statistical test according to them in the data file of theoretical resuts')
    parser.add_argument('--input_dir', '-i', help='directory that holds the results of he parametric bootstrapping', required=True)
    parser.add_argument('--theoretical_res_path', '-t', help='pah to the analysis results of the theoretical Tests', required=True)
    parser.add_argument('--output_path', '-o', help='path that will hold the results of the theoretical and empirical statistical Tests', required=True)
    parser.add_argument('--script_type', '-s', help='0 for RELAX, 1 for TraitRELAX', required=False, default=0)

    args = parser.parse_args()
    input_dir = args.input_dir
    theoretical_res_path = args.theoretical_res_path
    output_path = args.output_path
    script_type = int(args.script_type)

    integrated_df = pd.read_csv(theoretical_res_path)
    integrated_df['significant_by_PB'] = np.zeros(integrated_df.shape[0])
    integrated_df['pvalue_by_PB'] = np.zeros(integrated_df.shape[0])


    for path in os.listdir(input_dir):
        
        dataset_name = path

        print("beggining anaysis on dataset " + dataset_name)

        results_dir = input_dir + path + "/output/"

        if len(os.listdir(results_dir)) < 100:
            integrated_df.loc[integrated_df.dataset_id == dataset_name, 'significant_by_PB'] = "NA"
            continue

        if script_type == 0:
            script_path = "/groups/itay_mayrose/halabikeren/myScripts/python/bpp/extractRELAXResults.py"
        else:
            script_path = "/groups/itay_mayrose/halabikeren/myScripts/python/bpp/extractTraitRELAXResult.py"

        print("executing analysis of PB results on dataset " + dataset_name)
        res = os.system("python " + script_path + " -i " + results_dir + " -o " + results_dir + " > " + results_dir + "/mapping_job_to_dataset.log")

        print("integrating analysis result of PB on dataset " + dataset_name)

        results = pd.read_csv(results_dir + "res.csv")
        LR_scores = list(results["LRT_statistic"])
        LR_scores.sort() # sort in descending order
        LR_empirical_cutoff = LR_scores[int(len(LR_scores)*95/100)]
        integrated_df.loc[integrated_df.dataset_id == dataset_name, 'empirical_LRT_statistic_cutoff'] = LR_empirical_cutoff
        actual_LR = float(integrated_df.loc[integrated_df.dataset_id == dataset_name, 'LRT_statistic'])
        LR_scores.append(actual_LR)
        LR_scores.sort()
        pvalue_by_PB = 1 - LR_scores.index(actual_LR) / len(LR_scores)
        if actual_LR > float(LR_empirical_cutoff):
            integrated_df.loc[integrated_df.dataset_id == dataset_name, 'significant_by_PB'] = 1
        integrated_df.loc[integrated_df.dataset_id == dataset_name, 'pvalue_by_PB'] = pvalue_by_PB

    integrated_df.to_csv(output_path)