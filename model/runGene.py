
import os
import argparse

########################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--replicates', default=1, type=int, help="Replicates to run per simulation")
parser.add_argument('-g', '--geneID', required=True, type=str, help="Number of parallel simulations")
parser.add_argument('-od', '--outputdir', default='./', help="Output directory for simulation files")
parser.add_argument('--coTrsc', default=True, action=argparse.BooleanOptionalAction)
parser.add_argument('-t', '--simTime', default=60, type=int, help="Total biological time per replicate")
parser.add_argument('-wi', '--writeInterval', default=1.0, type=float, help="Biological time interval between frames written to trajectory file")
parser.add_argument('-gd', '--genomeDir', type=str)
parser.add_argument('-pd', '--paramDir', default=None)
parser.add_argument('-ps', '--paramSet', default=None)
parser.add_argument('-pf', '--paramFile', default=None)
parser.add_argument('-of', '--outputFile', default=None)
parser.add_argument('-k', '--kinetics', default='global')

args = parser.parse_args()
########################################################################

import lm
import isosplicer as splice
import hook

if not os.path.isdir(args.outputdir):
    os.mkdir(args.outputdir)

geneID = args.geneID

if (args.paramDir is not None) and (args.paramSet is not None):
    varsFile = f'{args.paramDir}/parameters_{args.paramSet}.xlsx'
elif args.paramFile is not None:
    varsFile = str(args.paramFile)
else:
    varsFile = None

simVars = splice.setSimVars(varsFile = varsFile)

sim = splice.initSim(simVars)
coTrsc = args.coTrsc

genes = {}

splice.getGeneData(args.genomeDir, geneID, genes)

splice.getIsoforms(sim, genes, args.genomeDir, geneID)

for geneID, geneDict in genes.items():
    
    splice.trscRxns(sim, genes, simVars, coTrsc=coTrsc)
    # splice.splicingRxns(sim, genes, simVars)
    if args.kinetics == 'global':
        splice.splicingRxnsUniform(sim, genes, simVars)
    elif args.kinetics == 'site':
        splice.splicingRxnsIndividualSites(sim, genes, simVars, varsFile)
    else:
        raise Exception(f"Options for kinetics are global or site, you provided: {args.kinetics}")


# if len(genes[geneID]['isoforms']) == 1:
#     replicates = 1
#     print('Single isoform, reducing replicates to 1')
# else:
replicates = args.replicates

if coTrsc:
    simTime = int(genes[geneID]['length']/int(simVars['Values']['trsc_elongation'])*60+args.simTime)
else:
    simTime = args.simTime

print(f'Simulation Time: {simTime}')
print(f'Write Interval: {args.writeInterval}')

sim.setWriteInterval(args.writeInterval)
sim.setSimulationTime(simTime)
sim.setHookInterval(1.0)

if args.outputFile is not None:
    filename = str(args.outputFile)
elif args.paramDir is not None:
    filename = f"{args.outputdir}{geneID}_{args.paramSet}.lm"
else:
    filename = f"{args.outputdir}{geneID}.lm"

os.system("rm -rf %s"%(filename))
sim.save(filename)

splicingSolver = hook.SplicingSolver(sim.species_id, genes, coTrsc = coTrsc)

# sim.runSolver(filename=filename, solver=splicingSolver, replicates=replicates, cudaDevices=None)

for rep in range(1,replicates+1):
    print(f'Replicate: {rep}')
    splicingSolver = hook.SplicingSolver(sim.species_id, genes, coTrsc = coTrsc, gparams=simVars)
    lm.runSolver(filename, solver=splicingSolver, replicate=rep, cudaDevices=[1], checkpointInterval=0)
    print()
