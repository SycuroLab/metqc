#!/usr/bin/python
import os
import sys
import re
import csv
import argparse
import subprocess

# Sample command
# python parse_snakemake_command_logs.py --log_infile /work/sycuro_lab/kevin/snakemake_pipeline_logs/metqc-analysis.log.txt --output_dir ~/testing_script/

parser = argparse.ArgumentParser()

log_infile = None
output_dir = None

parser.add_argument('--log_infile', action='store', dest='log_infile',
                    help='snakemake command log file as input. (i.e. *-analysis.log.txt)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME/output_dir)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

log_infile = results.log_infile
output_dir = results.output_dir

if(log_infile == None):
	print('\n')
	print('error: please use the --log_infile option to specify the snakemake command log file as input')
	print('log_infile =' + ' ' + str(log_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

# Snakemake pipeline success boolean flag. "true" for success and "false" for failure. Starts at "false" until script finds successful message.
snakemake_pipeline_success = "false"

# Get the filename of the log_infile.
log_basename = os.path.basename(log_infile)
#print(log_basename)

# Remove the ".log.txt" extension.
log_filename = log_basename.replace(".log.txt", "")
#print(log_filename)

#sys.exit()

# Dictionary of jobs where the key is the unique name consisting of the rule, sample id, and the slurm job id.
job_dict = {}

# List of rule components so that we can grab all the data from each line based on regex.
rule_list = []

# Variables for checking if snakemake finished without terminating.
(step_num_index, step_num_total, success_message, failure_message) = ("","","","")

# Open the log input file.
with open(log_infile, "r") as log_input_file:
    for line in log_input_file:
        line = line.rstrip('\n')
       
        # Parse the total number of steps.
        if(re.search("total             \d+", line)):
            rule_list = []
            match = re.search("total             (\d+)", line)# total number of steps.
            step_num_total = match.group(1)
            #print(line)
            #print(step_num_total)
            #sys.exit()
            #print(line)
            #print(rule_name)

        # Parse the name of the rule.
        if(re.search("rule (.*):$", line)):
            rule_list = []
            match = re.search("rule (.*):$", line)#rule metaphlan:
            rule_name = match.group(1)
            rule_list.append(rule_name)
            #print(line)
            #print(rule_name)

        # Parse the input parameters.
        if(re.search("    input: (.*)$", line)):
            match = re.search("    input: (.*)$", line)#    jobid: 34
            input_params = match.group(1)
            rule_list.append(input_params)

            #print(line)
            #print(input_params)

        # Parse the output parameters.
        if(re.search("    output: (.*)$", line)):
            match = re.search("    output: (.*)$", line)#    jobid: 34
            output_params = match.group(1)
            rule_list.append(output_params)

            #print(line)
            #print(output_params)

        # Parse the snakemake job id.
        if(re.search("    jobid: (\d+)$", line)):
            match = re.search("    jobid: (\d+)$", line)#    jobid: 34
            snakemake_job_id = match.group(1)
            rule_list.append(snakemake_job_id)
            
            #print(line)
            #print(snakemake_job_id)

        # Parse the sample id from the wildcards entry.
        if(re.search("wildcards: (sample=.*)$", line)):
            match = re.search("wildcards: (sample=.*)$", line)
            sample_id = match.group(1)
            rule_list.append(sample_id)
            
            #print(line)
            #print(sample_id)
        
                    
        # Submitted job 28 with external jobid 'Submitted batch job 10307744'.
        if(re.search("Submitted job (\d+) with external jobid \'Submitted batch job (\d+)\'\.$", line)):
            match = re.search("Submitted job (\d+) with external jobid \'Submitted batch job (\d+)\'\.$", line)
            snakemake_job_id = match.group(1)
            slurm_job_id = match.group(2)
            #print(snakemake_job_id)
            #print(slurm_job_id)
            
            # Execute the seff command on the slurm_job_id to get the resource usage for the job.
            result = subprocess.run(['seff', slurm_job_id], stdout=subprocess.PIPE)
            results_output = result.stdout.decode("utf-8")
            results_dict = {}
            for results_line in results_output.splitlines():
                #print(results_line)
                (key, value) = results_line.split(": ")
                results_dict[key] = value
            #print(results_dict)
            
            rule_list.append(results_dict)
            rule_list.append(snakemake_job_id)
            rule_list.append(slurm_job_id)
           
           
            (rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id) = ("","","","","","","","")
            
            # If the rule list is equal to 8 we have the following.
            if(len(rule_list) == 8):
                            
                (rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id) = rule_list
                #print(rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id)
                
                if(snakemake_job_id1 == snakemake_job_id2):
                    #print(rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id)
                    job_key = "|".join([slurm_job_id,rule_name,sample_id])
                    job_dict[job_key] = [rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id]
                    #print(job_key)
                    #print(job_dict[job_key])
            if(len(rule_list) == 7):
                            
                (rule_name,input_params,output_params,snakemake_job_id1,results_dict,snakemake_job_id2,slurm_job_id) = rule_list
                sample_id = "all"
                
                if(snakemake_job_id1 == snakemake_job_id2):
                    job_key = "|".join([slurm_job_id,rule_name,sample_id])
                    job_dict[job_key] = [rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id]
                    #print(job_key)
                    #print(job_dict[job_key])
            
        if(re.search("\d+ of \d+ steps \(100%\) done$", line)):
            match = re.search("((\d+) of (\d+) steps \(100%\) done$)", line)
            success_message = match.group(1)
            step_num_index = match.group(2)
            step_num_total = match.group(3)
            job_key = "|".join(["all",rule_name,"all"])
            job_dict[job_key] = [rule_name,input_params,output_params,snakemake_job_id1,"N/A","N/A","N/A","N/A"]
            snakemake_pipeline_success = "true"
            #print(step_num_index, step_num_total)

        if(re.search("Exiting because a job execution failed. Look above for error message$", line)):
            match = re.search("(Exiting because a job execution failed.) Look above for error message$", line)
            failure_message = match.group(1)

# Open the file for writing.
outfile = os.path.join(output_dir, "_".join([log_filename,"snakemake_job_resource_usage.log.txt"]))
output_file_handle = open(outfile, "w+")

# Resource efficiency.
# (snakemake) [kevin.muirhead@fc47 snakemake_pipeline_logs]$ seff 10248225
# Job ID: 10248225
# Cluster: arc-admin
# User/Group: kevin.muirhead/kevin.muirhead
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 7
# CPU Utilized: 00:08:53
# CPU Efficiency: 3.07% of 04:49:13 core-walltime
# Job Wall-clock time: 00:41:19
# Memory Utilized: 312.40 MB
# Memory Efficiency: 0.51% of 60.00 GB

# Number of steps counter.
num_steps_counter = 1

output_file_handle.write("\t".join(["Slurm Job ID", "Sample ID", "Rule Name", "Snakemake Job ID", "Cluster", "User/Group", "State", "Nodes", "Cores per node", "CPU Utilized", "CPU Efficiency","Job Wall-clock time","Memory Utilized","Memory Efficiency", "Input", "Output", "Step Number", "Percent Complete"]) + "\n")
for job_key in sorted(job_dict):
    rule_list = job_dict[job_key]
    print(rule_list)
    (rule_name,input_params,output_params,snakemake_job_id1,sample_id,results_dict,snakemake_job_id2,slurm_job_id) = rule_list
    if(results_dict != "N/A"):
        cluster_name = results_dict["Cluster"]
        user_and_group = results_dict["User/Group"]
        state = results_dict["State"]
        nodes = results_dict["Nodes"]
        cores_per_node = results_dict["Cores per node"]
        cpu_utilized = results_dict["CPU Utilized"]
        cpu_efficiency = results_dict["CPU Efficiency"]
        job_wall_time = results_dict["Job Wall-clock time"]
        memory_utilized = results_dict["Memory Utilized"]
        memory_efficiency = results_dict["Memory Efficiency"] 
    
    percent_complete = (float((num_steps_counter)/float(step_num_total)) * 100)
    #print(percent_complete)

    output_file_handle.write("\t".join([slurm_job_id, sample_id, rule_name, snakemake_job_id1, cluster_name, user_and_group, state, nodes, cores_per_node, cpu_utilized, cpu_efficiency, job_wall_time, memory_utilized, memory_efficiency, input_params, output_params, " ".join([str(num_steps_counter), "of", str(step_num_total)]), "{:0.2f}".format(percent_complete)]) + "\n")
    num_steps_counter += 1

    # Clear for the all job. The all job is not submitted as a cluster job.
    (cluster_name, user_and_group, state, nodes, cores_per_node, cpu_utilized, cpu_efficiency, job_wall_time, memory_utilized, memory_efficiency) = ("N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A","N/A")
output_file_handle.close()

if(snakemake_pipeline_success == "true"):
    print(log_filename + " " + "snakemake pipeline finished and completed successfully!!!")
    print(success_message)
else:
    print(log_filename + " " + "snakemake pipeline terminated and did not complete all the steps of the pipeline. Diagnose further.")
    print(failure_message)
    print("Check the following file for possible errors.")
    print(log_infile)
    
