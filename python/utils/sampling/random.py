import re, random

# auxiliary function that randomly chooses a subset from a sequence
def random_sample(input_path, ref_seq_name, subset_size, output_path):
    # catch all the sequence names in the fasta file
    sequences_catcher = re.compile(">([^\n]+)")
    with open(input_path, "r") as input:
        seqs = sequences_catcher.findall(input.read())
    sequences = []
    # remove the ">" sign from the beginning of each sequence name
    for seq in seqs:
        sequences.append(seq.replace(">", "", 1).replace(" ", ""))
    chosen_sequences = random.sample(sequences, subset_size)
    if ref_seq_name not in chosen_sequences:
        chosen_sequences.pop(1)
        chosen_sequences.append(ref_seq_name)
    # write the chosen sequences into a samples output file
    with open(output_path, "w") as output:
        for sequence in chosen_sequences:
            output.write(sequence + "\n")
    return chosen_sequences

# auxiliary function to extract the subset from the sample file
def parse_random_file(path):
    with open(path, "r") as file:
        subset = [s.replace("\n", "") for s in file.readlines()]
    return subset
