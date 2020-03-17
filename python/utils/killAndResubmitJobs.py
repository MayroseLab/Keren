import os, re

job_path_regex = re.compile("Submit_arguments\s*=\s-p\s\d*\s(.*?)project", re.DOTALL | re.MULTILINE)
error_path_regex = re.compile("Error_Path\s*=\s.*?\:(.*?)exec_host", re.DOTALL | re.MULTILINE)

if __name__ == '__main__':

    # get job IDs
    res = os.system("qstat -u halabikeren | awk '{print $1}' > /groups/itay_mayrose/halabikeren/jobs.txt")
    with open("/groups/itay_mayrose/halabikeren/jobs.txt", "r") as infile:
        jobs = infile.readlines()

    for job in jobs:
        if not ".power8." in job:
            jobs.remove(job)

    # jobs.remove('power8.tau.ac.il:\n')
    # jobs.remove('Job\n')
    # jobs.remove('\n')
    # jobs.remove("Req'd\n")
    # jobs.remove('---------------\n')

    for i in range(len(jobs)):
        if ".power8.\n" in jobs[i]:
            jobs[i] = jobs[i].replace(".power8.\n", "")

    print(jobs)
    exit()

    with open("jobs_information.csv", "w") as csv:
        csv.write("job_path,output_path\n")
        for job in jobs:
            print("cmd: ", "qstat -f " + job + " > /groups/itay_mayrose/halabikeren/job_info.txt")
            res = os.system("qstat -f " + job + " > /groups/itay_mayrose/halabikeren/job_info.txt")
            with open("/groups/itay_mayrose/halabikeren/job_info.txt", "r") as infile:
                content = infile.read()

            if "queue = itaym_others" in content or "queue = itaym3" in content:
                continue

            job_path = job_path_regex.search(content).group(1)
            error_path = error_path_regex.search(content).group(1)
            csv.write(job_path + "," + error_path + "\n")

            # # kill job
            # res = os.system("qdel " + job)
            #
            # # delete output files
            # output_file = error_path + job + ".power8.tau.ac.il.OU"
            # res = os.system("rm -r " + output_file)
            #
            # # resubmit the job
            # res = os.system("qsub -q itaym3 " + job_path)



