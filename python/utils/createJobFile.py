#### script for generation of .sh files
#  -*- coding: utf-8 -*-
__author__ = 'Shiran'
import socket, os


# set jobs environment
def set_job_env(job_files_path = "/groups/itay_mayrose/halabikeren/jobs/", error_files_path = "/groups/itay_mayrose/halabikeren/error_files/"):
    if not os.path.exists(job_files_path):
        res = os.system("mkdir -p " + job_files_path)
    if not os.path.exists(error_files_path):
        res = os.system("mkdir -p " + error_files_path)
    return 0

# creates a standard job file with a single command, with multi-threading option
def create_job_file(job_name, commands, file_name, error_files_path, job_files_path, priority, threads_number, touch_file_path, forbidden_nodes="compute-0-246", allowed_nodes=None,
                    limit_nodes=True, python=True, openmpi=True, language="tcsh", queue="itaymr", mem_alloc="2gb"):
    hostname = socket.gethostname()
    with open(job_files_path + "/" + file_name, "w") as handle:
        if language == "tcsh":
            handle.write("#!/bin/tcsh\n\n")  # Vladi said it should be tcsh!
            handle.write("#$ -N " + job_name + "\n")
            handle.write("#$ -S /bin/tcsh\n")
            handle.write("#$ -cwd\n")
            handle.write("#$ -l " + queue + "\n")
            handle.write("#$ -p " + str(priority) + "\n")
            handle.write("#$ -e " + error_files_path + "$JOB_NAME.$JOB_ID.ER\n")
            handle.write("#$ -o " + error_files_path + "$JOB_NAME.$JOB_ID.OU\n")
            handle.write("#$ -pe orte " + str(threads_number) + "\n")
            if allowed_nodes != None:
                handle.write("#$ -l h=" + allowed_nodes + "\n")
            elif limit_nodes:
                handle.write("#$ -l h=!" + forbidden_nodes + "\n")
            if openmpi:
                handle.write("module load rocks-openmpi\n")
            if python:
                if hostname == 'jekyl.tau.ac.il':
                    handle.write("module load python/python-3.3.0\n")
                elif hostname == 'lecs2.tau.ac.il':
                    handle.write("module load python/python-3.3.3\n")
                else:
                    handle.write("module load python/anaconda3-5.0.0\n")
        else:
            handle.write("# !/bin/bash\n\n")
            handle.write("#PBS -S /bin/bash\n")
            handle.write("#PBS -j oe\n")
            handle.write("#PBS -r y\n")
            handle.write("#PBS -q " + queue + "\n")
            handle.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
            handle.write("#PBS -N " + job_name + "\n")
            # handle.write("#PBS -k oe\n")
            handle.write("#PBS -e " + error_files_path + "\n")
            handle.write("#PBS -o " + error_files_path + "\n")
            if limit_nodes:
                handle.write("#PBS -l nodes=" + allowed_nodes[0] + "\n")
            else:
                handle.write("#PBS -l select=ncpus=1:mem=" + mem_alloc + "\n")
        for command in commands:
            if "python" not in command and "cd" not in command and threads_number > 1:
                handle.write("mpiexec -np " + str(threads_number) + " " + command + "\n")
            else:
                handle.write(command + "\n")
        handle.write("touch " + touch_file_path)
    return job_files_path + "/" + file_name


import os, re
def fix(indir):
    for file in os.listdir(indir):
        with open(indir+file, "r") as input:
            content = input.read()
        content = content.replace("model1.nodes_id", "true_history.model1.nodes_id")
        content = content.replace("model2.nodes_id", "true_history.model2.nodes_id")
        with open(indir+file, "w") as output:
            output.write(content)