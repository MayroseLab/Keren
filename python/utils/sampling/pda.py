import os, re, sys
from Bio import AlignIO, SeqIO
sys.path.insert(0, '/groups/itay_mayrose/halabikeren/python_scripts/')
from utils.alterTrees import get_tree_size

# hardcoded coefficient for the normalization value a to set the ratios of significance between PD score and weights score in PDA sampling
coeff = 1

# auxiliary function to extract the weighted pda PD's coefficient from the method's name
def get_wpda_coeff(method):
    coeff_catcher = re.compile("pda_(-*[0-9]+.?[0-9]*)")
    if "weight" not in method:
        return 1
    else:
        coeff_data = coeff_catcher.search(method)
    if coeff_data == None:
        return 1
    else:
        return float(coeff_data.group(1))

#### mapping functions ####

# auxiliary function that maps positions in a sequence to positions in a matching alignment of it
def align_seq_to_aln(sequence, alignment):
    seq_to_aln_map = dict()
    accumulated_gaps = 0
    for i in range(len(sequence)):
        j = i + accumulated_gaps
        while alignment[j] == '-':
            j += 1
            accumulated_gaps += 1
        seq_to_aln_map[i+1] = j+1
    return seq_to_aln_map

# auxiliary function that given two alignments with a common sequence, returns a mapping of one alignment to another according ot the common sequence
def map_ref_seq_positions(ref_seq_name, msa_path):
    # read the alignment of the reference sequence in msa
    sys.path.insert(0, '/groups/itay_mayrose/halabikeren/python_scripts/')
    from utils.alterAlignments import get_sequence_by_id
    reference_sequence = get_sequence_by_id(ref_seq_name, msa_path)
    reference_sequence_no_gaps = reference_sequence.replace("-", "")
    # map reference sequence positions to reference alignment positions according to msa
    ref_seq_to_ref_aln = align_seq_to_aln(reference_sequence_no_gaps, reference_sequence)
    return ref_seq_to_ref_aln

# auxiliary function to get th fraction of non-space characters in s given position over an alignment of sequences
def get_fraction(position, msa_path, chars=True):
    # print("get_fraction called with chars = " + str(chars)) # debug
    sys.path.insert(0, '/groups/itay_mayrose/halabikeren/python_scripts/')
    from utils.alterAlignments import get_format, get_alignment_sequences_amount
    msa_format = get_format(msa_path)
    seq_num = get_alignment_sequences_amount(msa_path)
    chars_counter = 0
    alignments = AlignIO.parse(open(msa_path), msa_format)
    alignment = next(alignments)
    for record in alignment:
        sequence = str(record.seq)
        if sequence[position-1] != "-":
            chars_counter += 1
    chars_fraction = float(chars_counter/seq_num)
    if not chars:
        return float(1-chars_fraction)
    return chars_fraction

# auxiliary function that gives taxon wright according to the informativity of its sequence in relations to the msa, using sum (i.e - favoring long sequences) or using avg (i.e - punishing long sequences)
def informativity_scoring(ref_seq_name, examined_seq_name, msa_path, position_to_coefficient):
    # fix sequence names formatting
    examined_seq_name = examined_seq_name.replace(" ", "")
    sys.path.insert(0, '/groups/itay_mayrose/halabikeren/python_scripts/')
    from utils.alterAlignments import get_sequence_by_id
    examined_sequence = get_sequence_by_id(examined_seq_name, msa_path)
    if len(examined_sequence) < 2:
        print("failed to extract examined sequence " + examined_seq_name)
        exit(1)
    examined_sequence_scores = []
    for position in list(position_to_coefficient.keys()):
        examined_sequence_pos_score = 0
        if examined_sequence[position-1] != "-":
            examined_sequence_pos_score = 1
        examined_sequence_scores.append(position_to_coefficient[position]*examined_sequence_pos_score)
    return examined_sequence_scores

# auxiliary function to return the mean of scoring vector
def informativity_avg_scoring(ref_seq_name, examined_seq_name, msa_path, position_to_coefficient, coeff=1):
    # fix sequence names formatting
    ref_seq_name = ref_seq_name.replace(" ", "")
    examined_seq_name = examined_seq_name.replace(" ", "")
    if ref_seq_name == examined_seq_name: # always favor the reference sequence
        return abs(coeff) * 10000000000
    informativity_scores = informativity_scoring(ref_seq_name, examined_seq_name, msa_path, position_to_coefficient)
    informativity_scores_avg = 0
    for score in informativity_scores:
        informativity_scores_avg += score
    return coeff * (informativity_scores_avg / len(informativity_scores))

# auxiliary function that gives default scoring to taxon as follows: it gives all taxa accept the ref taxon weight 1 and 42 to the ref taxon
def default_scoring(ref_seq_name, examined_seq_name, msa_path, position_to_coefficient, coeff=1):
    # fix sequence names formatting
    ref_seq_name = ref_seq_name.replace(" ", "")
    examined_seq_name = examined_seq_name.replace(" ", "")
    if ref_seq_name == examined_seq_name:
        return 42
    return 1

# auxiliary function that generates a file with taxon weights for PDA sampling procedure (in order to force PDA to choose the reference sequence)
def generate_weights_file(subset_size, tree_path, fasta_path, msa_path, ref_seq_name, weights_path, scoring_function=default_scoring, position_to_coefficient={}, coeff=1, debug=False):
    weights_map = dict()
    total_weights = 0
    seq_num = 0
    # remove spaces from req seq name
    ref_seq_name = ref_seq_name.replace(" ", "")
    # calculate the taxon's weights
    taxon_catcher = re.compile(">(.*)\n")
    with open(fasta_path, "r") as fasta:
        content = fasta.read()
        max_weight = -float('inf')
        for taxon in re.finditer(taxon_catcher, content):
            taxon_name = taxon.group(1)
            weights_map[taxon_name] = scoring_function(ref_seq_name, taxon_name, msa_path, position_to_coefficient, coeff=coeff)
            if taxon_name != ref_seq_name and max_weight < weights_map[taxon_name]:
                max_weight = weights_map[taxon_name]
            if taxon_name != ref_seq_name:
                total_weights += abs(weights_map[taxon_name])
        weights_map[ref_seq_name] = max_weight + abs(coeff)*len(list(position_to_coefficient.keys())) # make sure the ref seq is chosen by giving it maximal weight
        if debug:
            print("weights_map[ref_seq_name]: ", weights_map[ref_seq_name]) # debug
        total_weights += abs(weights_map[ref_seq_name])
    # calculate the rescaling parameter, in a manner that will give the plain PD score a proportional significance to the weights function, while considering the subsrt's size
    a_val = float(total_weights) / float(get_tree_size(tree_path)) # a_val is the coefficient of PD: WPD(G)=a_val*PD(G)+F(G)
    if debug:
        print("total_weights: ", float(total_weights))
        print("tree_size: ", float(get_tree_size(tree_path)))
        print("normalization_factor = ", a_val)
        print("weights_map: ", weights_map)
    # write the info into the weights file
    with open(weights_path, "w") as weights_file:
        # write the first line in the file: each branch receives rescaling size 1
        weights_file.write(str(a_val) + "\n")
        # write line for each taxon with its weight. all taxa but taxon 1 will receive weight 1, and taxon 1 will receive weight 42
        for taxon in list(weights_map.keys()):
            weights_file.write(taxon + " " + str(weights_map[taxon]) + "\n")
    return 0

# auxiliary function to run PDA
def run_pda(fasta_path, msa_path, tree_path, ref_seq_name, weights_path, output_path, subset_size, weights_method="default", coeff=1, debug=False, wpda_method=True):
    print("run_pda called with wpda_method value: " + str(wpda_method)) # debug
    # set the coring function according ot the weights method
    # get positions coefficients scores map
    position_to_coefficient = dict()
    if "weight" in weights_method:
        ref_seq_name = ref_seq_name.replace(" ", "")
        if ref_seq_name != "NA": # if there is a reference sequence, give the scoring based on the positions in the reference sequence only
            ref_seq_to_aln_map = map_ref_seq_positions(ref_seq_name, msa_path)
            if debug:
                print("ref_seq_to_aln_map: ", ref_seq_to_aln_map) # debug
            relevant_positions = list(ref_seq_to_aln_map.values()) # coefficients are calculated only on positions relevant to the reference sequence
            if debug:
                print("relevant_positions: ", relevant_positions) # debug
        else:  # else, give the scoring based on all the positions in the alignment
            sys.path.insert(0, "/groups/itay_mayrose/halabikeren/python_scripts/")
            from utils.alterAlignments import get_alignment_length
            relevant_positions = [i+1 for i in range(get_alignment_length(msa_path))]
        # get positions coefficients scores only on the positions of the reference sequence
        for position in relevant_positions:
            position_to_coefficient[position] = get_fraction(position, msa_path, chars=wpda_method)
        if debug:
            print("position_to_coefficient: ", position_to_coefficient) # debug
        scoring_function = informativity_avg_scoring
    else:
        scoring_function = default_scoring
    # generate a weights file according to the chosen scoring system
    res = generate_weights_file(subset_size, tree_path, fasta_path, msa_path, ref_seq_name, weights_path, scoring_function=scoring_function, position_to_coefficient=position_to_coefficient, coeff=coeff, debug=debug)
    # run PDA with the generated weights file
    res = os.system("pda " + tree_path + " " + output_path + " -k " + str(subset_size) + " -e " + weights_path)
    res = os.system("rm -r " + tree_path + ".log")
    return 0

# auxiliary function to extract from a pda file the sampled sequences' names
def parse_pda_file(pda_path):
    chosen_subset = []
    with open(pda_path, "r") as output:
        line = output.readline()
        while "The optimal PD" not in line: # skip all lines until the subset specification is reached
            line = output.readline()
        line = output.readline()
        while line != "" and line != "\n":              # add each sequence name to the chosen subset
            chosen_subset.append(line.replace("\n", ""))
            line = output.readline()
    return chosen_subset
