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

modes = ("critical", "single")
if len(argv) == 2 and argv[1] in modes:
    mode = argv[1]
else:
    raise ValueError("Unrecognised or missing calculation mode")


#dimensions = [1.5,1.7,1.8,1.9,2.,2.1,2.5,3.]
ns = [1.,1.2,1.4,1.6,1.8,2.]
options = [[{'-kappa': 1.5}],
           [{'-u': 0.075}],
           [{'-dim': 2}],
           [{'-rhomax': 3.5}],
           [{'-delta': 5e-4}],
           [{'-precision': 1e-6}],
           [{'-min_delta': 1e-4}],
           [{'-num_points': 200}],
           [{'-norm': 5}],
           [{'-a': 2}],
           [{'-N': n} for n in ns]]
           # [{'-d': 2.}, {'-d': 2.05}, {'-d': 2.07}, {'-d': 2.1}, {'-d': 2.2},
            # {'-d': 2.3}, {'-d': 2.5}, {'-d': 2.75}, {'-d': 2.9}, {'-d': 3}]]

all_configurations = options[0]
for single_option in options[1:]:
    new_configurations = []
    for single_config in single_option:
        for configuration_set in all_configurations:
            conf = single_config.copy()
            conf.update(configuration_set)
            if '-N' in conf.keys() and '-a' in conf.keys():
                conf['-out'] = 'flow_N=%1.1f_a=%i.csv' % (conf['-N'],conf['-a'])
            new_configurations.append(conf)
        
    all_configurations = new_configurations
project_dir = "/home/2/ac357729/Documents/FunctionalIntegrator"

config_file = open('%s/scripts/run.config' % project_dir, 'w')

if mode == "single":
    for conf in all_configurations:
        filename = "FunctionalIntegrator_single.sh"
        file = open("%s/scripts/%s" % (project_dir, filename), 'w')
        options = ["%s %s" % (k,v) for (k,v) in conf.items()]  
        
        basic_opt = [mode]
        opt = basic_opt + options
        
        opt_string = ' '.join(opt)
        config_file.write('%s: %s\n' % (filename, opt_string))
        
        file.write("cd %s\n" % project_dir)
        file.write("./FunctionalIntegrator %s\n" % opt_string)
        file.close()
        
        time = 120
        command = 'qsub -l mem=1000mb,walltime=%i:00:00 %s/scripts/%s' % (time, project_dir, filename)
        print command
        system(command)

if mode == "critical":
    for i, options in enumerate(all_configurations):
        filename = "FunctionalIntegrator_critical_%02i.sh" % i
        file = open("%s/scripts/%s" % (project_dir, filename), 'w')   
        
        basic_opt = [mode]
        opt = basic_opt + ['%s %s' % (k, v) for k, v in options.items()]
        
        opt_string = ' '.join(opt)
        
        config_file.write('%s: %s\n' % (filename, opt_string))
        
        file.write("cd %s\n" % project_dir)
        file.write("./FunctionalIntegrator %s\n" % opt_string)
        file.close()
    
        time = 120
        command = 'qsub -l mem=1000mb,walltime=%i:00:00 %s/scripts/%s' % (time, project_dir, filename)
        print command
        system(command)

