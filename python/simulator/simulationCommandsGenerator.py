import re, os

simulations_global_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/simulations/simulationStudyAug2019/"
jobs_dir = "/groups/itay_mayrose/halabikeren/jobs/TraitRELAX/simulations/simulationStudyAug2019/"
jobs_output_dir = "/groups/itay_mayrose/halabikeren/jobs_output/TraitRELAX/simulations/simulationStudyAug2019/"
relax_true_commands_path = simulations_global_dir + "relax_true_cmds.txt"
relax_mp_commands_path = simulations_global_dir + "relax_mp_cmds.txt"
traitrelax_commands_path = simulations_global_dir + "traitrelax_cmds.txt"
queue = "itaym1"
priority = 0

tbl_regex = re.compile("tbl_(\d*)_")
mu_regex = re.compile("mu_(\d*)_")
pi0_regex = re.compile("pi0_(\d*\.\d*)_")
kappa_regex = re.compile("kappa_(\d*\.?\d*)_")
p_regex = re.compile("p_(\d*\.?\d*)_")
omega1_regex = re.compile("omega1_(\d*\.?\d*)_")
omega2_regex = re.compile("omega2_(\d*\.?\d*)_")
theta1_regex = re.compile("theta1_(\d*\.?\d*)_")
theta2_regex = re.compile("theta2_(\d*\.?\d*)")

taxa_num_regex = re.compile("(\d*)_taxa")
pos_num_regex = re.compile("(\d*)_codons")
k_regex = re.compile("k_(\d*\.?\d*)")

if __name__ == '__main__':

    relax_true_jobs_dir = jobs_dir + "RELAXWithTrueHistory/"
    relax_true_jobs_output_dir = jobs_output_dir + "RELAXWithTrueHistory/"
    relax_mp_jobs_dir = jobs_dir + "RELAXWithMPHistory/"
    relax_mp_jobs_output_dir = jobs_output_dir + "RELAXWithMPHistory/"
    traitrelax_jobs_dir = jobs_dir + "TraitRELAX/"
    traitrelax_jobs_output_dir = jobs_output_dir + "TraitRELAX/"

    try:
        res=os.system("mkdir -p " + relax_true_jobs_dir)
        res = os.system("mkdir -p " + relax_true_jobs_output_dir)
        res=os.system("mkdir -p " + relax_mp_jobs_dir)
        res = os.system("mkdir -p " + relax_mp_jobs_output_dir)
        res=os.system("mkdir -p " + traitrelax_jobs_dir)
        res = os.system("mkdir -p " + traitrelax_jobs_output_dir)
    except:
        pass

    for path in os.listdir(simulations_global_dir):

        if not "tbl_" in path:
            continue

        tbl = int(tbl_regex.search(path).group(1))
        mu = float(mu_regex.search(path).group(1))
        pi0 = float(pi0_regex.search(path).group(1))
        kappa = float(kappa_regex.search(path).group(1))
        p = float(p_regex.search(path).group(1))
        omega1 = float(omega1_regex.search(path).group(1))
        omega0 = p * omega1
        omega2 = float(omega2_regex.search(path).group(1))
        theta1 = float(theta1_regex.search(path).group(1))
        p0 = theta1
        theta2 = float(theta2_regex.search(path).group(1))
        p1 = (1-theta1) * theta2

        indir = simulations_global_dir + path

        for inpath in os.listdir(indir):
            if not "taxa" in inpath:
                continue
            taxa_num = int(taxa_num_regex.search(inpath).group(1))
            inindir = indir + "/" + inpath + "/"

            # simulate trees
            trees_dir = inindir + "trees/"
            res = os.system('python /groups/itay_mayrose/halabikeren/myScripts/python/simulator/simulateTrees.py -o ' + trees_dir + ' -tn ' + str(taxa_num) + ' -mr ' + str(mu) + ' -ts ' + str(tbl) + ' -n 200')

            for ininpath in os.listdir(inindir):

                if "trees" in ininpath:
                    continue

                if not "codons" in ininpath:
                    continue

                pos_num = int(pos_num_regex.search(ininpath).group(1))
                ininindir = inindir + ininpath + "/"
                for inininpath in os.listdir(ininindir):
                    k = float(k_regex.search(inininpath).group(1))
                    replicates_num = 50
                    executions_num = 50
                    if k == 1:
                        executions_num = 200
                        replicates_num = 200

                    # simulate data
                    simulations_dir = ininindir + inininpath + "/"
                    cmd='python /groups/itay_mayrose/halabikeren/myScripts/python/simulator/TraitRELAXSimulator.py -t ' + trees_dir + ' -mu ' + str(mu) + ' -pi0 ' + str(pi0) + ' -kappa ' + str(kappa) + ' -omega0 ' + str(omega0) + ' -omega1 ' + str(omega1) + ' -omega2 ' + str(omega2) + ' -p0 ' + str(p0) + ' -p1 ' + str(p1) + ' -rep ' + str(replicates_num) + ' -al ' + str(pos_num) + ' -k ' + str(k) + ' -o ' + simulations_dir
                    print(cmd)
                    res=os.system(cmd)

                    # write maximum parsimony histories
                    relax_true_parameters_dir = simulations_dir + "relax_param/"
                    relax_mp_parameters_dir = simulations_dir + "mp_param/"
                    traitrelax_parameters_dir = simulations_dir + "traitrelax_param/"
                    cmd = 'python /groups/itay_mayrose/halabikeren/myScripts/python/bpp/CreateMPHistories.py -t ' + trees_dir + ' -c ' + simulations_dir + ' -p ' + relax_true_parameters_dir + ' -o ' + relax_mp_parameters_dir
                    print(cmd)
                    res=os.system(cmd)

                    # create jobs and jobs output directories
                    relax_true_jobs_subdir = relax_true_jobs_dir + ininindir + "/" + inininpath
                    if not os.path.exists(relax_true_jobs_subdir):
                        res = os.system("mkdir -p " + relax_true_jobs_subdir)
                    relax_true_jobs_output_subdir = relax_true_jobs_output_dir + ininindir + "/" + inininpath
                    if not os.path.exists(relax_true_jobs_output_subdir):
                        res = os.system("mkdir -p " + relax_true_jobs_output_subdir)
                    relax_mp_jobs_subdir = relax_mp_jobs_dir + ininindir + "/" + inininpath
                    if not os.path.exists(relax_mp_jobs_subdir):
                        res = os.system("mkdir -p " + relax_mp_jobs_subdir)
                    relax_mp_jobs_output_subdir = relax_mp_jobs_output_dir + ininindir + "/" + inininpath
                    if not os.path.exists(relax_mp_jobs_output_subdir):
                        res = os.system("mkdir -p " + relax_mp_jobs_output_subdir)
                    traitrelax_jobs_subdir = traitrelax_jobs_dir + ininindir + "/" + inininpath
                    if not os.path.exists(traitrelax_jobs_subdir):
                        res = os.system("mkdir -p " + traitrelax_jobs_subdir)
                    traitrelax_jobs_output_subdir = traitrelax_jobs_output_dir + ininindir + "/" + inininpath
                    if not os.path.exists(traitrelax_jobs_output_subdir):
                        res = os.system("mkdir -p " + traitrelax_jobs_output_subdir)

                    # now write the relevant commands into the 3 command files
                    relax_true_cmd = 'python "/groups/itay_mayrose/halabikeren/myScripts/python/bpp/runBppProgram.py" -pp "/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/relax" -pd ' + relax_true_parameters_dir + ' -jd ' + relax_true_jobs_subdir + ' -err ' + relax_true_jobs_output_subdir + ' -rn ' + str(executions_num) + ' -q ' + queue + ' -pr ' + str(priority)
                    with open(relax_true_commands_path, "a+") as infile:
                        infile.write(relax_true_cmd)

                    relax_mp_cmd = 'python "/groups/itay_mayrose/halabikeren/myScripts/python/bpp/runBppProgram.py" -pp "/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/relax" -pd ' + relax_mp_parameters_dir + ' -jd ' + relax_mp_jobs_subdir + ' -err ' + relax_mp_jobs_output_subdir + ' -rn ' + str(executions_num) + ' -q ' + queue + ' -pr ' + str(priority)
                    with open(relax_mp_commands_path, "a+") as infile:
                        infile.write(relax_mp_cmd)

                    traitrelax_cmd = 'python "/groups/itay_mayrose/halabikeren/myScripts/python/bpp/runBppProgram.py" -pp "/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/traitrelax" -pd ' + traitrelax_parameters_dir + ' -jd ' + traitrelax_jobs_subdir + ' -err ' + traitrelax_jobs_output_subdir + ' -rn ' + str(executions_num) + ' -q ' + queue + ' -pr ' + str(priority)
                    with open(traitrelax_commands_path, "a+") as infile:
                        infile.write(traitrelax_cmd)







