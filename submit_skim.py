#!/usr/binA/env python                                                                                                                 
# IMPORT MODULES NEEDED                                                                                                               
import argparse
import glob
import logging
import os
import subprocess
import time
import shlex

# How to run example:
# python3 submit_skim.py -i  /home/vamitamas/Samples8GeV/v3.3.3_ecalPN-batch3/ -py selnskim.py -s "sim"

def main():
    #READ IN AND PROCESS ARGUMENTS
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--indir', action='store', dest='indir')
    parser.add_argument('-t', '--test', action='store_true', dest='test')
    parser.add_argument('-py', '--script', action='store', dest='script')
    parser.add_argument('-s', '--suffix', action='store', dest='suffix')
    args = parser.parse_args()

    indirs = glob.glob(args.indir + "/*.root")

    logging.basicConfig(format='[ submitJobs ][ %(levelname)s ]: %(message)s', level=logging.DEBUG)


    #MOVE TO INDIR, CREATE JOB AND LOG DIRECTORY
    os.chdir(args.indir)
    jobdir = '/sdf/home/t/tamasvami/jobs/'
    #logdir = os.getcwd()+'/logs/'
    logdir = '/sdf/home/t/tamasvami/logs/'
    if not os.path.exists(jobdir):
        logging.info('Creating jobs directory.')
        os.makedirs(jobdir)
    if not os.path.exists(logdir):
        logging.info('Creating log files directory.')
        os.makedirs(logdir)

    batch_command = 'sbatch -p milano'
    command = 'singularity run --no-home --bind /sdf/home/t/tamasvami:/sdf/home/t/tamasvami:rw,/fs/ddn/sdf/group/ldmx:/fs/ddn/sdf/group/ldmx:rw %s . fire %s ' % (os.environ.get('LDMX_PRODUCTION_IMG'),args.script) #not including file name yet


    #WRITE JOB FILES AND SUBMIT
    jobs = len(indirs)
    for job in range(0,jobs):
        specific_command = command + '%s' % (indirs[job])
        if args.suffix:
            specific_command += " -s " + args.suffix
        with open('%s/slurm_submit_%d.job' % (jobdir,job), 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('#SBATCH --nodes=1 --ntasks-per-node=1\n')
            f.write('#SBATCH --error=%s/slurm-%%A_%%a.err\n' % logdir)
            f.write('#SBATCH --output=%s/slurm-%%A_%%a.out\n' % logdir)
            f.write('#SBATCH --time=00:03:00\n\n')                                                      
            #f.write('#SBATCH --mail-type=END,FAIL\n')
            #f.write('#SBATCH --mail-user=tamasvami@ucsb.edu\n')
            f.write('cd $SLURM_SUBMIT_DIR\n\n')
            f.write('/bin/hostname\n\n')
            f.write('%s\n' % specific_command)

        submit_command = '%s %s/slurm_submit_%d.job' % (batch_command, jobdir, job)
        logging.info('Job submit command: %s' % submit_command)
        reform_specific = specific_command.replace("/", "\/")
        if args.test:
            if job==jobs-1:
                logging.info('Test command: %s' % specific_command)
                #subprocess.Popen(specific_command, shell=True).wait()
                #subprocess.run("cp /home/aminali/sizeskim.py .", shell=True)
                #print(reform_specific)
                #subprocess.Popen(reform_specific, shell=True)
                #print(["singularity run"]+specific_command.split()[2:])
                #subprocess.run(["singularity run"]+specific_command.split()[2:], shell=True)
                #subprocess.run("ls \/home\/aminali\/sizeskim.py", shell=True)
        else:
            subprocess.Popen(submit_command, shell=True).wait()
            time.sleep(2)

if __name__ == "__main__" :
    main()

