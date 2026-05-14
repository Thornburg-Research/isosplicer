import threading
import subprocess
from subprocess import Popen, PIPE

import concurrent.futures

import os

import argparse

########################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--replicates', default=1, type=int, help="Replicates to run per simulation")
parser.add_argument('-th', '--threads', default=1, type=int, help="Number of parallel simulations")
parser.add_argument('-od', '--outputdir', default='./', help="Output directory for simulation files")
parser.add_argument('--coTrsc', default=True, action=argparse.BooleanOptionalAction)
parser.add_argument('-t', '--simTime', default=60, type=int, help="Total biological time per replicate")
parser.add_argument('-wi', '--writeInterval', default=1.0, type=float, help="Biological time interval between frames written to trajectory file")
parser.add_argument('-gl', '--geneListFile', type=str, help="File containing list of genes to be simulated")
parser.add_argument('-gd', '--genomeDir', type=str)
parser.add_argument('-pd', '--paramDir', default=None)
parser.add_argument('-ps', '--parameterSets', default=None)
parser.add_argument('-k', '--kinetics', default='global')

args = parser.parse_args()
########################################################################


########################################################################
def runSim(args, geneID, taskID, paramSet=None):

    if args.coTrsc:
        sim_args = ['python', 'runGene.py', '-r', str(int(args.replicates)), '-od', args.outputdir, '--coTrsc', '-g', geneID, '-dd', args.genomeDir, '-t', str(int(args.simTime)), '-wi', str(float(args.writeInterval)), '-k', str(args.kinetics)]
        
    elif not args.coTrsc:
        sim_args = ['python', 'runGene.py', '-r', str(int(args.replicates)), '-od', args.outputdir, '--no-coTrsc', '-g', geneID, '-dd', args.genomeDir, '-t', str(int(args.simTime)), '-wi', str(float(args.writeInterval)), '-k', str(args.kinetics)]

    if paramSet is not None:
        if args.paramDir is None:
            raise ValueError("Directory for parameter files must be provided when running multiple parameter sets.")
        sim_args.append('-pd')
        sim_args.append(args.paramDir)
        sim_args.append('-ps')
        sim_args.append(str(int(paramSet)))

    print(sim_args)
    if not args.coTrsc:
        print('No Cotranscriptional Splicing')

    if paramSet is None:
        log_file = f'{args.outputdir}logs/log_{geneID}.log'
    else:
        log_file = f'{args.outputdir}logs/log_{geneID}_{paramSet}.log'

    with open(log_file, "w") as log:

        subprocess.run(sim_args, stdout=log, stderr=subprocess.STDOUT, text=True)
    # ,
    #                stdout=log_file,
    #                stderr=log_file,
    #                text=True,  # To handle text-based output
    #                check=True,  # Raise an exception for non-zero exit code
    #               )

    print(f'Completed splicing for gene: {geneID}')

    return taskID
########################################################################


########################################################################
if __name__ == "__main__":

    if not os.path.isdir(args.outputdir):
        os.mkdir(args.outputdir)
    os.mkdir(f'{args.outputdir}logs')

    with open(args.geneListFile, 'r') as f:
        geneList = [line.rstrip() for line in f]

    total_tasks = int(len(geneList))

    max_workers = int(args.threads)

    print(f'Running {total_tasks} simulations with maximum concurrent simulations {max_workers}')

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit initial batch of tasks

        if args.parameterSets is None:
            futures = {executor.submit(runSim, args, geneList[i], i, None): i for i in range(total_tasks)}
        else:
            total_tasks = total_tasks*int(args.parameterSets)
            print(f'Total tasks: {total_tasks}')
            futures = {executor.submit(runSim, args, geneList[int(i//int(args.parameterSets))], i, int(i%int(args.parameterSets)+1)): i for i in range(total_tasks)}
        
        # futures = {executor.submit(runSim, args, geneList[i], i): i for i in range(total_tasks)}

        for future in concurrent.futures.as_completed(futures):
            
            taskID = futures[future]
            
            try:
                result = future.result()
                if args.parameterSets is None:
                    print(f"Simulation {result}: Gene {geneList[result]} completed successfully")
                else:
                    print(f"Simulation {result}: Gene {geneList[int(result//int(args.parameterSets))]} parameter set {int(result%int(args.parameterSets)+1)} completed successfully")
                    
            except Exception as exc:
                print(f"Gene {taskID} generated an exception: {exc}")

########################################################################


