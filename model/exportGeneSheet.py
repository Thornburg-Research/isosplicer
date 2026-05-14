
import os
import argparse

########################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--geneID', required=True, type=str, help="Number of parallel simulations")
parser.add_argument('-od', '--outputdir', default='./', help="Output directory for simulation files")
parser.add_argument('--coTrsc', default=True, action=argparse.BooleanOptionalAction)
parser.add_argument('-gd', '--genomeDir', type=str)
parser.add_argument('-exp', '--export', default='global')
parser.add_argument('-pb', '--proteinBinders', default=0, type=int)

args = parser.parse_args()
########################################################################

import lm
import isosplicer as splice

if not os.path.isdir(args.outputdir):
    os.mkdir(args.outputdir)

geneID = args.geneID

simVars = splice.setSimVars()

sim = splice.initSim(simVars)
coTrsc = args.coTrsc

genes = {}

splice.getGeneData(args.genomeDir, geneID, genes)

splice.getIsoforms(sim, genes, args.genomeDir, geneID)

if args.export == 'global':
    splice.exportGlobalParameters(genes, simVars, outputDir=args.outputdir)
elif args.export == 'site':
    splice.exportSiteParameters(genes, simVars, outputDir=args.outputdir, proteinBinders=args.proteinBinders)
