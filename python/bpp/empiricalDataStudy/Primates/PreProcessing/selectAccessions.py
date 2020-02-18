import os, argparse, pickle, pandas as pd
from Bio import Entrez, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML


def map_prot_to_nuc_accession(species_to_accessions):
    Entrez.email = 'halabikeren@gmail.com'
    prot_accessions_nuc_accessions = dict()
    for species in species_to_accessions:
        for accession in species_to_accessions[species]:
            request = Entrez.epost("nucleotide", id=accession)
            result = Entrez.read(request)
            webEnv = result["WebEnv"]
            queryKey = result["QueryKey"]
            handle = Entrez.efetch(db="nucleotide", retmode="xml", webenv=webEnv, query_key=queryKey)
            for r in Entrez.parse(handle):
                for item in r['GBSeq_feature-table']:
                    relevant_items = []
                    if "GBFeature_quals" in item and "protein_id" in str(item) and ("SEMG2" in str(item) or "semenogelin II" in str(item) or "semenogelin 2" in str(item)):
                        relevant_items.append(item["GBFeature_quals"])
                    for item in relevant_items:
                        for element in item:
                            if "GBQualifier_name" in element and element["GBQualifier_name"] == "protein_id":
                                prot_accessions_nuc_accessions[element["GBQualifier_value"]] = accession

    return prot_accessions_nuc_accessions

def get_accessions_evaluations(BLAST_result_path, prot_accessions_nuc_accessions):
    accession_to_evaluation = dict()
    results = open(BLAST_result_path, 'r')
    records = NCBIXML.parse(results)
    for record in records:
        for alignment in record.alignments:
            blast_accession = (alignment.accession).lstrip()
            try:
                nuc_accession = prot_accessions_nuc_accessions[blast_accession + ".1"]
                hsp = alignment.hsps[0]
                score = hsp.score
                expect = hsp.expect
                # identities = hsp.identities
                # positives = hsp.positives
                gaps = hsp.gaps
                align_length = hsp.align_length
                coverage = (align_length - gaps) / align_length
                accession_to_evaluation[nuc_accession] = (score, expect, coverage)
            except Exception as e:
                pass
    return accession_to_evaluation

def get_best_accession(accessions_list, accession_to_evaluation, nuc_to_prot_accessions):
    # assign any accession that blast returned as hit as initial accession
    best_accession = accessions_list[0]
    for accession in accessions_list:
        if accession in accession_to_evaluation and accession in nuc_to_prot_accessions:
            best_accession = accession
        break
    if best_accession == "":  # if no accession of the spece was returned by blast, return NA
        return "NA"

    # check if there is any beter accession, and if there is, return it
    for accession in accessions_list:
        try:
            if accession in accession_to_evaluation:
                score = accession_to_evaluation[accession][0]
                evalue = accession_to_evaluation[accession][1]
                coverage = accession_to_evaluation[accession][2]
                if score > accession_to_evaluation[best_accession][0]:
                    best_accession = accession
                elif score == accession_to_evaluation[best_accession][0]:
                    if coverage > accession_to_evaluation[best_accession][2]:
                        best_accession = accession
                    elif coverage == accession_to_evaluation[best_accession][2]:
                        if evalue < accession_to_evaluation[best_accession][1]:
                            best_accession = accession
        except:
            pass

    return best_accession

def _get_best_verified_accession(accessions_list, accession_to_evaluation, prot_accessions_nuc_accessions):

    # get the protein accessions of each accession in the list, if available
    # print("accessions_list: ", accessions_list)
    nuc_to_prot_accessions = dict()
    for nuc_accession in accessions_list:
        if nuc_accession in prot_accessions_nuc_accessions.values():
            for prot_accession in prot_accessions_nuc_accessions:
                if prot_accessions_nuc_accessions[prot_accession] == nuc_accession:
                    nuc_to_prot_accessions[nuc_accession] = prot_accession
    # print("nuc_to_prot_accessions: ", nuc_to_prot_accessions)

    best_accession = get_best_accession(accessions_list, accession_to_evaluation, nuc_to_prot_accessions)

    # BLASTP NR while using the best protein accession as query and make sure the majority of the protein accessons of the other nucleotide accessions in the list are returned as hits
    best_verified = False
    while not best_verified and len(accessions_list) > 1:
        matched_accessions = []
        fraction_of_matched_accessions = 0
        try:
            result_handle = NCBIWWW.qblast("blastp", "refseq_protein", nuc_to_prot_accessions[best_accession])
        except:
            accessions_list.remove(best_accession)
            if len(accessions_list) > 0:
                best_accession = get_best_accession(accessions_list, accession_to_evaluation, nuc_to_prot_accessions)
            else:
                best_verified = True
            continue
        records = NCBIXML.parse(result_handle)

        for record in records:
            for alignment in record.alignments:
                matched_accessions.append(alignment.accession)
        for nuc_accession in nuc_to_prot_accessions:
            prot_accession = (nuc_to_prot_accessions[nuc_accession]).replace(".1","")
            if prot_accession in matched_accessions:
                fraction_of_matched_accessions += 1
        fraction_of_matched_accessions = fraction_of_matched_accessions / len(nuc_to_prot_accessions.keys())
        if fraction_of_matched_accessions < 0.5:
            print("best accession by score failed verification")
            accessions_list.remove(best_accession)
            if len(accessions_list) > 0:
                best_accession = get_best_accession(accessions_list, accession_to_evaluation, nuc_to_prot_accessions)
            else:
                best_verified = True
        else:
            best_verified = True
    return best_accession

def run_mafft(input_path, output_path):
    res = os.system("module load module load mafft/mafft-7.407")
    res = os.system("mafft --localpair --maxiterate 1000 " + input_path + " > " + output_path)
    if res != 0:
        return -1
    return 0

def get_similarity_score(sequence_1, sequence_2):
    similarity = 0
    num_of_pos = len(sequence_1)
    for pos in range(num_of_pos):
        if sequence_1[pos] == sequence_2[pos]:
            similarity += 1
    similarity = similarity / num_of_pos
    return similarity

def get_best_verified_accession(accessions_list, prot_accessions_nuc_accessions, species, msas_dir):

    # if there is only one accession in the list - chose it
    if len(accessions_list) == 1:
        return accessions_list[0]

    # MSA the protein accessions and take the one with the highest similarity to the others
    sequences_dir = msas_dir + species.replace(" ", "_") + "/"
    if not os.path.exists(sequences_dir):
        res = os.system("mkdir -p " + sequences_dir)

    # get the protein accessions of each accession in the list, if available
    nuc_to_prot_accessions = dict()
    for nuc_accession in accessions_list:
        if nuc_accession in prot_accessions_nuc_accessions.values():
            for prot_accession in prot_accessions_nuc_accessions:
                if prot_accessions_nuc_accessions[prot_accession] == nuc_accession:
                    nuc_to_prot_accessions[nuc_accession] = prot_accession

    # create a file with the protein sequences of all the accessions in the list
    with open(sequences_dir + "unaligned.fas", "w") as infile:
        for nuc_accession in accessions_list:
            try:
                prot_accession = nuc_to_prot_accessions[nuc_accession]
                Entrez.email = "halabikeren@gmail.com"
                request = Entrez.epost("protein", id=prot_accession)
                result = Entrez.read(request)
                webEnv = result["WebEnv"]
                queryKey = result["QueryKey"]
                handle = Entrez.efetch(db="protein",retmode="xml", webenv=webEnv, query_key=queryKey)
                records = [r for r in Entrez.parse(handle)]
                prot_sequence = records[0]["GBSeq_sequence"].upper()
                infile.write(">" + nuc_accession + "\n" + prot_sequence + "\n") # name the sequence by its nucleotide accession because this is what we are looking for
            except Exception as e:
                # print("error processing protein sequence: ", e)
                continue


    # align the sequences with mafft
    run_mafft(sequences_dir + "unaligned.fas", sequences_dir + "aligned.fas")

    # measure the total similarity of each accession against the rest and report the one with the highest similarity
    alignment = AlignIO.read(sequences_dir + "aligned.fas", "fasta")
    accession_to_similarity = dict()
    records = [record for record in alignment]
    for record in records:
        name = record.id
        total_similarity = 0
        for other_record in records:
            if record.id != other_record.id:
                total_similarity += get_similarity_score(record.seq, other_record.seq)
        accession_to_similarity[name] = total_similarity

    best_accession = list(accession_to_similarity.keys())[0]
    for accession in accession_to_similarity:
        if accession_to_similarity[accession] > accession_to_similarity[best_accession]:
            best_accession = accession

    return best_accession







if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Creates a union between two sequence files based on the species names in the sequences identifications')
    parser.add_argument('--BLAST_result_path', '-i', help='full path for xml of BLAST result', required=False,
                        default="/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/primate_SEMG2/sequenceDataProcessing/PSI-BLASTP-Result.xml")
    parser.add_argument('--msas_dir', '-m', help='directory to hold the MSA of competing accessions for each specie', required=False, default="/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/primate_SEMG2/sequenceDataProcessing/competingAccessionsMSAs/")
    parser.add_argument('--output_path', '-o', help='path to a list that matches to each specie an accession', required=False, default="/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/primate_SEMG2/sequenceDataProcessing/competingAccessionsMSAs/chosen_accessions.csv")
    parser.add_argument('--sequence_data_output_path', '-s', help='path to hold the unaligned sequences', required=False, defalut="/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/primate_SEMG2/sequenceDataProcessing/unalignedSeqData.fas")

    args = parser.parse_args()
    BLAST_result_path = args.BLAST_result_path
    msas_dir = args.msas_dir
    output_path = args.output_path
    sequence_data_output_path = args.sequence_data_output_path

    species_to_accessions = {'Pan troglodytes':['AY259287.1', 'AY781386.1', 'DP000037'],
                            'Gorilla gorilla':['AY259289.1', 'AY781387.1', 'DP000041'],
                            'Pan paniscus':['AY259288.1'],
                            'Pongo pygmaeus':['AY259290.1', 'AY781388.1', 'DP000045'],
                            'Pongo abelii':['NM_001168574.1'],
                            'Symphalangus syndactylus':['LC148360.1', 'LC148359.1', 'LC148358.1', 'LC148357.1', 'LC148356.1', 'LC148355.1', 'LC148354.1', 'LC148353.1', 'LC148352.1', 'LC148351.1', 'LC148350.1'],
                            'Nomascus leucogenys':['LC148349.1', 'LC148348.1', 'LC148347.1'],
                            'Hylobates lar':['AY781389.1', 'LC148339.1',  'LC148338.1', 'LC148337.1', 'LC148336.1', 'LC148335.1', 'LC148334.1', 'LC148333.1', 'LC148332.1', 'LC148331.1', 'LC148330.1', 'LC148329.1', 'LC148328.1', 'LC148327.1', 'LC148326.1', 'LC148325.1', 'LC148324.1','LC148323.1', 'LC148322.1', 'LC148321.1', 'LC148320.1', 'LC148319.1', 'LC148318.1', 'LC148317', 'LC148316.1', 'LC148315.1', 'LC148314.1', 'LC148313.1', 'LC148312.1','LC148311.1'],
                            'Nomascus gabriellae':['LC148346.1'],
                            'Macaca nemestrina':['AY781391.1'],
                            'Colobus guereza':['AY781392.1', 'DP000038'],
                            'Nasalis larvatus':['LC216902.1'],
                            'Hylobates pileatus':['LC148345.1', 'LC148344.1', 'LC148343.1', 'LC148342.1', 'LC148341.1', 'LC148340.1'],
                            'Mandrillus sphinx':['LC216901.1'],
                            'Hylobates agilis':['LC148310.1', 'LC148309.1'],
                            'Papio anubis':['DP000036'],
                            'Callithrix jacchus':['DP000044'],
                            'Lemur catta':['DP000042'],
                            'Ateles geoffroyi':['AY781393.1'],
                            'Homo sapiens':['AY259284.1', 'AY259285.1', 'AY259286.1', 'NM_003008.3'],
                            'Macaca fascicularis':['AY781390.1', 'AY781390.1'],
                            'Aotus nancymaae':['DP000046'],
                            'Macaca mulatta':['DP000043'],
                            'Chlorocebus aethiops':['DP000048','X92589.1'],
                            'Hylobates klossii':['AY259291.1'],
                            'Otolemur_garnettii':['DP000040']}

    # translate the nucleotide annotations into protein annotations
    if not os.path.exists(msas_dir + "prot_accessions_nuc_accessions.txt"):
        prot_accessions_nuc_accessions = map_prot_to_nuc_accession(species_to_accessions)
        with open(msas_dir + "prot_accessions_nuc_accessions.txt", "wb") as outfile:
            pickle.dump(prot_accessions_nuc_accessions, outfile, protocol=4)
    else:
        with open(msas_dir + "prot_accessions_nuc_accessions.txt", "rb") as infile:
            prot_accessions_nuc_accessions = pickle.load(infile)
    # print("prot_accessions_nuc_accessions: ", prot_accessions_nuc_accessions)

    # process the blast result: for each accession, save the following: accessin (key), values (from the liagnments hsp): score, expect, align_length  (length of alignment), gaps (number of gaps), positives (#aligned positions)
    # there are two records from PSI-BLAST - for me, only the first record is relevant (if an accession doesn't appear in it - it's too different from the homo spaiens sequence and should be excluded)
    # in each record there are alignments. Each alignment instance represents one hit, for which the statistics are saved in the hsps data memeber in the first index. the accession is saved in the field accession
    # BLAST_result_path = "C:/Users/ItayMNB7/Desktop/PSI-BLASTP-Result.xml"
    # accession_to_evaluation = get_accessions_evaluations(BLAST_result_path, prot_accessions_nuc_accessions)
    # print("# relevant hits = ", len(accession_to_evaluation))

    # for each specie, report the accession with the best score
    data = []
    for species in species_to_accessions:
        best_accession = get_best_verified_accession(species_to_accessions[species], prot_accessions_nuc_accessions, species, msas_dir)
        data.append({"species": species, "chosen_accession": best_accession})
    df = pd.DataFrame(data)
    df.to_csv(output_path)

    # write the sequence data into output file
    with open(sequence_data_output_path, "w") as outfile:
        for species in species_to_accessions:
            chosen_accession = list(df.loc[df['species'] == species]['chosen_accession'])[0]
            Entrez.email = "halabikeren@gmail.com"
            request = Entrez.epost("nucleotide", id=chosen_accession)
            result = Entrez.read(request)
            webEnv = result["WebEnv"]
            queryKey = result["QueryKey"]
            handle = Entrez.efetch(db="nucleotide", retmode="xml", webenv=webEnv, query_key=queryKey)
            for r in Entrez.parse(handle):
                for item in r['GBSeq_feature-table']:
                    relevant_items = []
                    if "GBFeature_quals" in item and "CDS" in str(item) and ("SEMG2" in str(item) or "semenogelin II" in str(item) or "semenogelin 2" in str(item)):
                        relevant_items.append(item["GBFeature_quals"])
                    for item in relevant_items:
                        for element in item:
                            if "GBQualifier_name" in element and element["GBQualifier_name"] == "protein_id":
                                prot_accessions_nuc_accessions[element["GBQualifier_value"]] = accession

            outfile.write(">" + species + "\n" + sequence + "\n")


