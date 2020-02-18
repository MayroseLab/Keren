import sys
sys.path.insert(0, '/groups/itay_mayrose/halabikeren/python_scripts/utils/')
from alterAlignments import get_alignment_sequences_amount

# auxiliary function to intialize a dictionary with subsets sampled by each combination of fraction and method
def initialize_samples_map(sampling_fractions, sampling_methods):
    samples_map = dict()
    for fraction in sampling_fractions:
        samples_map[fraction] = dict()
        for method in sampling_methods:
            samples_map[fraction][method] = "NA"
    return samples_map

# auxiliary function to initialize the required subset sizes
def initialize_samples_sizes(fasta_path, sampling_fractions):
    dataset_size = get_alignment_sequences_amount(fasta_path)
    samples_sizes = dict()
    for fraction in sampling_fractions:
        samples_sizes[fraction] = int(dataset_size * fraction / 100)
    return samples_sizes

# auxiliary function to initialize a map of the redcued inputs paths
def initialize_reduced_inputs_map(sampling_fractions, sampling_methods):
    reduced_input_paths = dict()
    for fraction in sampling_fractions:
        reduced_input_paths[fraction] = dict()
        for method in sampling_methods:
            reduced_input_paths[fraction][method] = dict()
            reduced_input_paths[fraction][method]["fasta"] = "NA"
            reduced_input_paths[fraction][method]["msa"] = "NA"
            reduced_input_paths[fraction][method]["tree"] = "NA"
    return reduced_input_paths
