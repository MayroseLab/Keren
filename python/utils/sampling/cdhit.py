import os, re

# auxiliary function to set the given garbage directory as global
def set_garbage_dir(dir):
    if not os.path.exists(dir):
        res = os.system("mkdir " + dir)
    output_dir = dir

# Recursive function that given an interval of thresholds [top, bottom] amd their compatible k values, finds the compatible threshold for k in the interval
def find_threshold(input_file, k, top, top_k, bottom, bottom_k, k_dict, output_dict, garbage_path):

    if (top - bottom) < 0.01:
        print("there is no convergence of thresholds for the clusters number: " + str(k))
        print("the closest clusters' number to the required one is: " + str(top_k))
        return top_k

    if top <= bottom:  # stop condition, in case no compatible threshold for the required k was found
        return top_k
    if k == top_k:
        return top_k
    if k == bottom_k:
        return bottom_k

    internal = ((top - bottom) / ((top_k - bottom_k) / (k - bottom_k))) % 1 + bottom  # calculate the nternal threshold for the next iteration using the gap between the thresholds' interval and thei comatible k's intervalG
    internal_k, cdhit_output = count_clusters(input_file, internal, garbage_path)
    k_dict[internal_k] = internal
    output_dict[internal_k] = cdhit_output
    if k in k_dict:  # stop condition, in case the compatible threshold for the required k was found
        return k

    elif internal_k < k:  # if internal_k < k -> the internal threshold is not strict enough -> continue with the interval [internal, top]
        return find_threshold(input_file, k, top, top_k, internal, internal_k, k_dict, output_dict, garbage_path)
    else:  # if internal_k > k -> the internal threshold is too strict -> continue with the interval [bottom,internal]
        return find_threshold(input_file, k, internal, internal_k, bottom, bottom_k, k_dict, output_dict, garbage_path)

# Given a fasta file and a threshold between 0.4 and 1, the function returns the number of clusters given by running cd-hit on the fasta file with the given threshold
def count_clusters(input_file, threshold, garbage_path):
    output_file = run_cd_hit(input_file, threshold, garbage_path)
    clusters = open(output_file + '.clstr')
    clusters_number = 0
    for line in clusters:
        if line.startswith(">Cluster"):
            clusters_number = clusters_number + 1
    clusters.close()
    return clusters_number, output_file

# auxiliary function to run cd-hit
def run_cd_hit(input_file, threshold, output_dir):
    input_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = output_dir + input_name + '_threshold-' + str(threshold)
    if threshold > 0.7:
        wordlen = '5'
    elif threshold > 0.6:
        wordlen = '4'
    elif threshold > 0.5:
        wordlen = '3'
    else:
        wordlen = '2'
    os.system('cd-hit -i ' + input_file + ' -o ' + output_file + ' -c ' + str(threshold) + ' -n ' + wordlen)
    return output_file

# given a fasta file and a required clusters number, find the threshold that causes cd-hit to returns k clusters
def find_optimal_threshold(input_file, k, garbage_path):

    k_dict = dict()  # k_dict holds k as key and its compatible threshold as value
    output_dict = dict() # output_dict holds the names of the cd-hit files produced for each examined threshold

    top_k, cdhit_output = count_clusters(input_file, 1, garbage_path)  # calculate the k returned by cd-hit run with thershold=1 (top threshold)
    k_dict[top_k] = 1
    output_dict[top_k] = cdhit_output
    if k > top_k:
        print("there is no threshold big enough for this clusters number")
        return top_k, k_dict[top_k], cdhit_output

    bottom_k, cdhit_output = count_clusters(input_file, 0.4, garbage_path)  # calculate the k returned by cd-hit run with thershold=0.4 (bottom threshold)
    k_dict[bottom_k] = 0.4
    output_dict[bottom_k] = cdhit_output
    if k < bottom_k:
        print("there is no threshold small enough for this clusters number")
        print("closest clusters num is: " + str(bottom_k) + " with the threshold: " + str(k_dict[bottom_k]))
        return bottom_k, k_dict[bottom_k], cdhit_output

    res = find_threshold(input_file, k, 1, top_k, 0.4, bottom_k, k_dict, output_dict, garbage_path)  # calculate the threshold compatible for thr required k using the help function find_threshold
    return res, k_dict[res], output_dict[res]

# auxiliary function to extract the chosen subset from a cd-hit output file
def parse_cdhit_file(cdhit_path, ref_seq_name):
    chosen_subset = []
    seq_catcher = re.compile(">([^\n]+)")
    with open(cdhit_path, "r") as output:
        seqs = seq_catcher.findall(output.read())
    for seq in seqs:
        chosen_subset.append(seq.replace(">", "", 1).replace(" ", ""))
    # exchange the cluster representative of the ref sequence with the ref seqence itself
    if ref_seq_name != "NA":
        ref_representative = get_ref_cluster_representative(cdhit_path, ref_seq_name)
        chosen_subset.append(ref_seq_name)
        chosen_subset.remove(ref_representative)
    return chosen_subset

# auxiliary function that finds the cluster that holds the refrence sequence, and returns its representative
def get_ref_cluster_representative(input_file, ref_seq_name):
    if ref_seq_name != "NA":
        input_path = input_file + ".clstr"
        ref_cluster_catcher = re.compile(">Cluster\s([0-9]+)\n.*>" + ref_seq_name + "...")
        refs_cluster_cacther = re.compile(">Cluster\s([0-9]+)\n.*>" + ref_seq_name)
        representative_catcher = re.compile(">(.+)\.\.\. \*")
        with open(input_path, "r") as input:
            input_content = input.read()
            representatives = representative_catcher.findall(input_content)
            refs_cluster = refs_cluster_cacther.search(input_content).group(1)
            refs_representative = representatives[int(refs_cluster)]
        return refs_representative
    else:
        return None

