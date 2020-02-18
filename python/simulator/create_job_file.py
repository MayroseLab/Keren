__author__ = 'Shiran'
import socket


def create_indelible_job_file(job_name, output_dir, control_file, file_name, error_files_path, job_files_path):
    with open(job_files_path + "/" + file_name, "w") as handle:
        handle.write("#!/bin/tcsh\n\n")  # Vladi said it should be tcsh!
        handle.write("#$ -N " + job_name + "\n")
        handle.write("#$ -S /bin/tcsh\n")
        handle.write("#$ -cwd\n")
        handle.write("#$ -l itaym\n")
        handle.write("#$ -e " + error_files_path + "$JOB_NAME.$JOB_ID.ER\n")
        handle.write("#$ -o " + error_files_path + "$JOB_NAME.$JOB_ID.OU\n")
        handle.write("cd " + output_dir + "\n")
        handle.write("echo " + control_file + " | /groups/pupko/haim/Programs/indelible/INDELibleV1.03/src/indelible\n")
        handle.write("touch flag_indelible_done")
    return job_files_path + "/" + file_name


def create_job_file(job_name, command, file_name, error_files_path, job_files_path):
	hostname = socket.gethostname()
	with open(job_files_path + "/" + file_name, "w") as handle:
		handle.write("#!/bin/tcsh\n\n")  # Vladi said it should be tcsh!
		handle.write("#$ -N " + job_name + "\n")
		handle.write("#$ -S /bin/tcsh\n")
		handle.write("#$ -cwd\n")
		handle.write("#$ -l itaym\n")
		handle.write("#$ -e " + error_files_path + "$JOB_NAME.$JOB_ID.ER\n")
		handle.write("#$ -o " + error_files_path + "$JOB_NAME.$JOB_ID.OU\n")
		if hostname == 'jekyl.tau.ac.il':
			handle.write("module load python/python-3.3.0\n")
		elif hostname == 'lecs2.tau.ac.il':
			handle.write("module load python/python-3.3.3\n")
		handle.write(command + "\n")
	return job_files_path + "/" + file_name
