import os, pickle, argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import translate


# get the length of the sequences in the alignment
def get_alignment_length(input_path, input_format="fasta"):
    with open(input_path, "r") as input_file:
        alignment = AlignIO.read(input_file, input_format)
    return alignment.get_alignment_length()


# translate codon alignment to aa alignment
def translate_msa(codon_path, aa_path):
    input_handle = open(codon_path, "r")
    output_handle = open(aa_path, "w")
    alignment = SeqIO.parse(input_handle, "fasta")
    while True:
        try:
            record = next(alignment)
            seq_name = record.id
            codon_seq = record.seq
            gap_pos = []
            # print("seq: ", seq_name)
            # print("content: ", codon_seq)
            for i in range(len(codon_seq)//3):
                if codon_seq[3*i] == "-":
                    gap_pos.append(i)
            # print("gap_pos: ", gap_pos)
            aa_seq_content = []
            for i in range(len(codon_seq)//3):
                if i not in gap_pos:
                    aa_seq_content.append(str(translate(codon_seq[3*i:3*i+3])))
                else:
                    aa_seq_content.append("-")
            seq_content = "".join(aa_seq_content)
            output_handle.write(">" + seq_name + "\n")
            output_handle.write(seq_content + "\n")
        except Exception as e:
            print("error: ", e)
            break
    output_handle.close()
    input_handle.close()
    return 0


# auxiliary function to generate an old to new positions mapping
def get_pos_map(map_info_path, msa_path):
    num_of_pos = get_alignment_length(msa_path)
    # print("num_of_pos: ", num_of_pos)
    pos_map = dict()
    with open(map_info_path, "r") as map_info:
        map_info.readline()
        info = map_info.readline()
        saved_pos = [int(i)+1 for i in info.split(",")]
        # print("saved_pos: ", saved_pos)
    accumulated_gaps = 0
    for pos in range(1, num_of_pos+1): # fix 27.4.17 - positions here are counted as codon positions
        # print("pos: ", pos)
        if pos in saved_pos:
            pos_map[pos] = pos - accumulated_gaps
            # print("matched pos: ",  pos - accumulated_gaps)
        else:
            pos_map[pos] = 0
            # print("matched pos: 0")
            accumulated_gaps += 1
    return pos_map


# auxiliary function to reduce a single sequence by the indices of the non-informative positions
def reduce_seq(seq_content, relevant_pos, debug=False):
    seq_content = seq_content.strip()
    if debug:
        print("reduce_seq called")
        print("seq_content: ")
        print(seq_content)
    new_seq_cont = []
    for pos in relevant_pos:
        if pos == 1:
            base_pos = 0
        else:
            base_pos = 3*(pos-1)
        if debug:
            print("required range: (" + str(base_pos) + ", " + str(base_pos+3) + ")") # debug
            print("content in it: " + seq_content[base_pos:base_pos+3]) # debug
        new_seq_cont.append(seq_content[base_pos:base_pos+3])
    return "".join(new_seq_cont)


# auxiliary function that translates the codons in the alignment to amino acids
def reverse_translate_msa(codons_reference_path, pos_map, codon_output_path, format="fasta"):
    records = list(SeqIO.parse(codons_reference_path, format))
    new_pos = [i for i in pos_map.keys() if pos_map[i] != 0]
    with open(codon_output_path, "w") as outfile:
        for record in records:
            seq_name = record.id
            aa_seq = record.seq
            codon_seq = reduce_seq(str(aa_seq), new_pos, debug=False)
            outfile.write(">" + seq_name + "\n" + codon_seq + "\n")
    return 0


# function that filters out non-informative codon positions with trimal
# assumes input format is fasta
def filter_by_trimal(codon_msa_path, aa_msa_path, trimmed_aa_msa_path, trimmed_codon_msa_path, map_path, trimal_cutoff):
    
    # if the alignment is made of codons -> get a compatible msa with amino acids instead of codons
    orig_pos_num = get_alignment_length(aa_msa_path)
    semi_pos_map = map_path + "_semi_map_info" # generate AAs positions map
    res = os.system("/groups/itay_mayrose/halabikeren/programs/trimAl/source/trimal -in " + aa_msa_path + " -out " + trimmed_aa_msa_path + " -gt " + str(trimal_cutoff) + " -colnumbering > " + semi_pos_map)
    if res != 0:
        return -1

    pos_map = get_pos_map(semi_pos_map, aa_msa_path)
    with open(map_path, "wb") as map_file:
        pickle.dump(pos_map, map_file, protocol=0)
    res = reverse_translate_msa(codon_msa_path, pos_map, trimmed_codon_msa_path)
    trimmed_pos_num = get_alignment_length(trimmed_aa_msa_path)
    filtered_pos_fraction = (orig_pos_num - trimmed_pos_num) / orig_pos_num
    if filtered_pos_fraction > 0.3:
        print("Please note: when trimmed with cuttof " + str(trimal_cutoff) + " more than 30% of the positions have been filtered out")
    return 0

if __name__ == '__main__':
    
    # process input from command line
    parser = argparse.ArgumentParser(description='trims positions from a codon alignment based on the respective protein alignment')
    parser.add_argument('--original_codon_msa_path', '-ic', help='path to the input codon alignment', required=True)
    parser.add_argument('--original_aa_msa_path', '-ia', help='path to the input aa alignment', required=True)
    parser.add_argument('--trimmed_codon_msa_path', '-oc', help='path to the trimmed codon alignment', required=True)
    parser.add_argument('--trimmed_aa_msa_path', '-oa', help='path to the trimmed aa alignment', required=True)
    parser.add_argument('--trimmed_positions_map', '-om', help='path to the map of trimmed aa positions', required=True)
    parser.add_argument('--cutoff', '-c', help='the amount of gap positions required for an aa position to be filtered out', required=True)

    args = parser.parse_args()
    original_codon_msa_path = args.original_codon_msa_path
    original_aa_msa_path = args.original_aa_msa_path
    trimmed_codon_msa_path = args.trimmed_codon_msa_path
    trimmed_aa_msa_path = args.trimmed_aa_msa_path
    trimmed_positions_map = args.trimmed_positions_map
    cutoff = float(args.cutoff)

    translate_msa(original_codon_msa_path, original_aa_msa_path)
    filter_by_trimal(original_codon_msa_path, original_aa_msa_path, trimmed_aa_msa_path, trimmed_codon_msa_path, trimmed_positions_map, cutoff)
