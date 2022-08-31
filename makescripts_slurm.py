from os import system
from sys import argv

# def options_dict(dim):
#     at3D = 0.2
#     at2D = 1.9
#     if dim < 2:
#         k = at2D
#         r = 1.025
#     elif dim <= 3:
#         k = (3 - dim) * at2D + (dim - 2) * at3D
#         r = (3 - dim) * 1.025 + (dim - 2) * 10
#     else:
#         k = at3D
#         r = 10
#
#     return {'-dim': dim, '-out': 'flow_d_%1.2f.csv' % (dim), '-kappa': k, '-rhomax': r}

test = False

modes = ("critical", "single")
if len(argv) == 2 and argv[1] in modes:
    mode = argv[1]
else:
    raise ValueError("Unrecognised or missing calculation mode")

# options = [[{'-kappa': .2}],
#            [{'-u': .5}],
#            [{'-dim': 3}],
#            [{'-rhomax': 12.}],
#            [{'-delta': 1e-5}],
#            [{'-precision': 1e-6}],
#            [{'-num_points': 200}],
#            [{'-norm': 5}],
#            [{'-a': i * .1}  for i in range(10, 31)],
#            [{'-N': 2}],
#            [{'-sigma_normalization': 'false'}]]
#

options = [[{'-kappa': .86 + i * 0.005} for i in range(0, 12)],
           [{'-u': 0.075}],
           [{'-dim': 2}],
           [{'-rhomax': 1.6}],
           [{'-delta': 5e-6}],
           [{'-precision': 1e-6}],
           [{'-norm_point': 0}],
           [{'-num_points': 160}],
           [{'-a': 1.7 + i * .1}  for i in range(0, 5)],
           [{'-N': 2}],
           [{'-sigma_normalization': 'false'}]]

# options = [[{'-kappa': .50 + i * 0.02} for i in range(0, 4)],
#            [{'-u': 1}],
#            [{'-dim': 2}],
#            [{'-rhomax': 1.2}],
#            [{'-delta': 5e-6}],
#            [{'-precision': 1e-6}],
#            [{'-num_points': 120}],
#            [{'-a': i * .1}  for i in range(10, 22)],
#            [{'-N': 2}],
#            [{'-sigma_normalization': 'false'}]]

all_configurations = options[0]
for single_option in options[1:]:
    new_configurations = []
    for single_config in single_option:
        for configuration_set in all_configurations:
            conf = single_config.copy()
            conf.update(configuration_set)
            if '-N' in conf.keys():
                conf['-out'] = 'fixed_norm_flow_k=%f_a=%f.csv' % (conf['-kappa'], conf['-a'])
            new_configurations.append(conf)
        
    all_configurations = new_configurations

project_dir = "/clusteruy/home/achlebicki/FunctionalIntegrator"

config_file = open('%s/scripts/run.config' % project_dir, 'w')

MAX_RUNNING_JOBS = 30.
MAX_CPUS = 80.
MAX_TOTAL_JOBS = 50

def submit_jobs():
    if len(all_configurations) > MAX_TOTAL_JOBS:
        print "Number of submitted jobs is bigger than the limit"
        
        decision = None
        while decision not in ("", "n", "N", "y", "Y"):
            decision = raw_input("Do you want to continue (y/N) ")
        
        if decision in ("", "n", "N"):
            print "Job submission aborted"
            return
    
    
    for i, conf in enumerate(all_configurations):
        time = 60
        # if mode == "single":
        thread_count = int((i + 1) * MAX_CPUS / MAX_RUNNING_JOBS) - int(i * MAX_CPUS / MAX_RUNNING_JOBS)
        
        conf['-threads'] = thread_count
    
        opt_list = ["%s %s" % (k, v) for (k, v) in conf.items()]  
        
        options = ' '.join([mode] + opt_list)
        filename = "%s_%02i.sh" % (mode, i)
        file = open("%s/scripts/%s" % (project_dir, filename), 'w')   
        
        config_file.write('%s: %s\n' % (filename, options))
        
        file.write("cd %s\n" % project_dir)
        file.write("./FunctionalIntegrator %s\n" % options)
        file.close()
            
        filename = "%s_%02i.sh" % (mode, i)
        file = open("%s/scripts/%s" % (project_dir, filename), 'w')
        file.write("#!/bin/bash \n")
        file.write("#SBATCH --job-name=%s \n" % filename)
        file.write("#SBATCH --ntasks=1\n")
        file.write("#SBATCH --cpus-per-task=%i \n" % thread_count)
        file.write("#SBATCH --mem=128 \n")
        file.write("#SBATCH --time=%i:00:00 \n" % time)
        file.write("#SBATCH --tmp=1G \n")
        file.write("#SBATCH --partition=normal \n")
        file.write("#SBATCH --qos=normal \n")
        file.write("#SBATCH --mail-type=ALL \n")
        #file.write("#SBATCH --mail-user=achlebicki@fing.edu.uy \n\n\n")
        file.write("source /etc/profile.d/modules.sh \n\n\n")
        file.write("echo Number of nodes: $SLURM_JOB_NUM_NODES \n")
        file.write("echo Node list: $SLURM_NODELIST \n")
        file.write("echo Number of tasks: $SLURM_TASKS_PER_NODE \n")
        file.write("echo -ne 'Execution starts on ';date \n")
        file.write("cd %s \n" % project_dir)
        file.write("./FunctionalIntegrator %s\n\n" % options)
        file.write("echo -ne 'Execution ends ';date \n")
        
        config_file.write('%s: %s\n' % (filename, options))
        
        file.close()
    
        command = 'sbatch %s/scripts/%s' % (project_dir, filename)
        print(command)
        
        if not test:
            system(command)

submit_jobs()
