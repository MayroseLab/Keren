import os, re, pickle
from Bio import AlignIO, SeqIO
from Bio.Seq import translate


amino_acid_chars = ["A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N","E","F","I","L","P","Q"]

# translate the names of the original sequences (taxa) to plain number, and return the new msa path and the mapping of names to numbers
def convert_sequences_names(fasta_path, updated_fasta_path, prefix="S"):
    with open(fasta_path, "r") as fasta_file:
        content = fasta_file.read()
    new_content = content.replace(">", ">"+prefix)
    with open(updated_fasta_path, "w") as new_fasta_file:
        new_fasta_file.write(new_content)
    return 0

# creates a new msa file with the new names according ot the old to new names translator
def update_sequences_names(input_msa, output_msa, names_translator_path):
    with open(names_translator_path, "rb") as names_translator_file:  # load the names map
        names_translator = pickle.load(names_translator_file)
    with open(input_msa, "r") as in_msa:
        msa_content = in_msa.read()
    for old_name in list(names_translator.keys()):
        msa_content = msa_content.replace(old_name, names_translator[old_name], 1)
    with open(output_msa, "w") as out_msa:
        out_msa.write(msa_content)
    return 0

# get the type of characters in the simulation (amino acids or nucleotides)
def get_msa_type(msa_path, input_format):
    alignment = AlignIO.read(open(msa_path), input_format)
    record = alignment[0]
    content = record.seq
    for aa in amino_acid_chars:
        if aa in content:
            return "AA"
    return "NUC"

# get the number of sequences in an alignment, depending on its format
def get_alignment_sequences_amount(input_path):
    input_format = get_format(input_path)
    if input_format == "fasta":
        alignment = SeqIO.parse(open(input_path), input_format)
    else:
        alignment = AlignIO.parse(open(input_path), input_format)
    seq_num = 0
    while True:
        try:
            record = next(alignment)
            seq_num += 1
        except:
            return seq_num

# get the length of the sequences in the alignment
def get_alignment_length(input_path):
    input_format = get_format(input_path)
    alignment = AlignIO.read(open(input_path), input_format)
    return alignment.get_alignment_length()

# generate mafft lindsi alignment: input must be in fasta format, output is given in fasta format
def run_mafft(input_path, output_path):
    if not os.path.exists(output_path):
        res = os.system("mafft --localpair --maxiterate 1000 " + input_path + " > " + output_path)
        if res != 0:
            return -1
    return 0

# nomi was here
def get_format(path):
    with open(path, 'r') as f:
        line = f.readline()
    if ">" in line:
        return "fasta"
    return "phylip"

# auxiliary function that receives input fasta, output path and a subset, and copies the sequences in the suvbset from the input fasta to the output path
def reduce_file(input_path, output_path, subset):
    format = get_format(input_path)
    alignment = list(SeqIO.parse(open(input_path), format))
    reduced_alignment = [record for record in alignment if record.id in subset]
    with open(output_path, "w") as output_file:
        SeqIO.write(reduced_alignment, output_file, format)
    return 0

# auxiliary function that extracts an alignment of a sequence by its id from a given alignment
def get_sequence_by_id(sequence_id, msa_path):
    msa_format = get_format(msa_path)
    with open(msa_path, "r") as msa:
        alignments = SeqIO.parse(msa, msa_format)
        while True:
            try:
                record = next(alignments)
                if str(record.id) == sequence_id:
                    return str(record.seq)
            except:
                return ""
    return ""


#### convert from [fasta, phylip, nexus] to [fasta, phylip, nexus] ####

def convert_fasta_to_phylip(input_path, output_path, blank=False):
    if blank:
        intermediata_path = input_path.replace(".fas", ".translated_fasta")
        names_translator_path = input_path.replace(".fas", ".names_map")
        res = convert_sequences_names(input_path, intermediata_path, names_translator_path, src="FastML")
        input_handle = open(intermediata_path, "rU")
    else:
        input_handle = open(input_path, "rU")
    output_handle = open(output_path, "w")
    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "phylip-relaxed")
    input_handle.close()
    output_handle.close()
    return 0

def convert_phylip_to_fasta(input_path, output_path):
    input_handle = open(input_path, "rU")
    output_handle = open(output_path, "w")
    alignments = AlignIO.parse(input_handle, "phylip-relaxed")
    AlignIO.write(alignments, output_handle, "fasta")
    output_handle.close()
    input_handle.close()
    return 0

def convert_fasta_to_nexus(input_path, output_path):
    res = os.system("perl /groups/pupko/elilevy/spartaABC/PS/PS_pipe/scripts/fasta2nexus.pl " + input_path + " " + output_path)
    return 0

def convert_phylip_to_newick(input_path, output_path):
    middle = input_path.replace("phy", "aln")
    convert_phylip_to_fasta(input_path, middle)
    convert_fasta_to_nexus(middle, output_path)
    return 0

def convert_formats(input_path, output_path, input_format, output_format):
    input_handle = open(input_path, "rU")
    output_handle = open(output_path, "w")
    alignments = AlignIO.parse(input_handle, input_format)
    AlignIO.write(alignments, output_handle, output_format)
    output_handle.close()
    input_handle.close()
    return 0

#### filter non-informative codon positions with trimal ####

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
def translate_aas_to_codons(codons_reference_path, pos_map, codon_output_path):
    # print("pos_map: ", pos_map)
    # print("codons_reference_path: ", codons_reference_path)
    # print("codon_output_path: ", codon_output_path)
    seq_cather = re.compile("(>[0-9]+)([^>]+)")
    with open(codons_reference_path, "r") as codons_ref:
        codon_ref_content = codons_ref.read()
    codon_output_file = open(codon_output_path, "w")
    new_pos = [i for i in pos_map.keys() if pos_map[i] != 0]
    # print("new_pos: ", new_pos)
    for seq in seq_cather.finditer(codon_ref_content):
        seq_name = seq.group(1).replace(">", ">S")
        old_seq_content = seq.group(2)
        # print("old_seq_content: ", old_seq_content)
        new_seq_content = reduce_seq(old_seq_content.replace("\n", ""), new_pos, debug=False)
        # print("new_seq_content: ", new_seq_content)
        # print("")
        codon_output_file.write(seq_name + "\n")
        codon_output_file.write(new_seq_content + "\n")
    codon_output_file.close()
    return 0

def translate_codons_to_aas(codon_path, aa_path):
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
        except:
            break
    output_handle.close()
    input_handle.close()
    return 0

# function that filters out non-informative codon positions with trimal
def filter_by_trimal(msa_path, output_path, map_path, trimal_cutoff, codons=True):
    # if the msa is in phylip format, convert it to fasta format
    if ".phy" in msa_path:
        codons_fasta_msa_path = msa_path.replace("phy", "fas")
        res = convert_phylip_to_fasta(msa_path, codons_fasta_msa_path)
    # if the alignment is made of codons -> get a compatible msa with amino acids instead of codons
    else:
        codons_fasta_msa_path = msa_path + ".copy"
        res = os.system("cp " + msa_path + " " + codons_fasta_msa_path)
    # print("codons_fasta_msa_path: ", codons_fasta_msa_path)
    aas_fasta_msa_path = codons_fasta_msa_path.replace("fas", "aas_fas")
    # print("aas_fasta_msa_path: ", aas_fasta_msa_path)
    aas_fasta_output_path = output_path.replace("fas", "aas_fas")
    # print("aas_fasta_output_path: ", aas_fasta_output_path)
    res = translate_codons_to_aas(msa_path, aas_fasta_msa_path)
    orig_pos_num = get_alignment_length(aas_fasta_msa_path)
    semi_pos_map = map_path + "_semi_map_info" # generate AAs positions map
    res = os.system("/groups/pupko/haim/Programs/trimAl/source/trimal -in " + aas_fasta_msa_path + " -out " + aas_fasta_output_path + " -gt " + str(trimal_cutoff) + " -colnumbering > " + semi_pos_map)
    if res != 0:
        return -1
    # update the output file to have correct names
    trimmed_pos_num = get_alignment_length(aas_fasta_msa_path)
    updated_aas_fasta_path = aas_fasta_output_path.replace(".aas_fas", "_updated.aas_fas")
    updated_fasta_output_file = open(updated_aas_fasta_path, "w")
    name_catcher = re.compile("(>[0-9]+)[^\n]+")
    # res = os.rename(aas_fasta_output_path, updated_aas_fasta_path) # rename didn't suffice here because trimal adds stuff to the sequences names
    with open(aas_fasta_output_path, "r") as in_file:
        for line in in_file:
            if ">" in line:
                full_name = name_catcher.search(line)
                updated_fasta_output_file.write(full_name.group(1)+"\n")
            else:
                updated_fasta_output_file.write(line)
        updated_fasta_output_file.close()
    updated_codons_fasta_path = aas_fasta_output_path.replace("aas_fas", "fas")
    pos_map = get_pos_map(semi_pos_map, aas_fasta_msa_path)
    with open(map_path, "wb") as map_file:
        pickle.dump(pos_map, map_file, protocol=0)
    res = translate_aas_to_codons(msa_path, pos_map, updated_codons_fasta_path)
    # res = convert_fasta_to_phylip(updated_codons_fasta_path, output_path)
    filtered_pos_fraction = (orig_pos_num - trimmed_pos_num) / orig_pos_num
    if filtered_pos_fraction > 0.3:
        print("Please note: when trimmed with cuttof " + trimal_cutoff + " more than 30% of the positions have been filtered out")
    res = os.system("rm -r " + codons_fasta_msa_path)
    res = os.system("rm -r " + aas_fasta_output_path)
    res = os.system("rm -r " + updated_aas_fasta_path)
    res = os.system("rm -r " + semi_pos_map)
    return 0


