# isosplicer

## Genome Data

Before running the Python executables here, you need to downlaod and extract the human genome information with the instructions in ```human_genome/```

This extraction only needs to be done once, it is not required for each new installation of isosplicer.

## Running Instructions

### Environment

You need to be in the ```isosplicer``` conda environment described in the parent directory to run the Python executables in this directory.

```
conda activate isosplicer
```

### ```exportGeneSheet.py```: Exporting an editable parameters file

| Argument | Short | Default | Description |
|----------|------|---------|-------------|
| `--geneID` | `-g` | - | Name of gene from human genome |
| `--genomeDir` | `-gd` | - | Directory containing extracted human genome files (e.g. SRSF6.xlsx) |
| `--outputdir` | `-od` | `"./"` | Directory to save output Excel file |
| `--export` | `-exp` | `"global"` | Export global (global) or global+site-specific (site) parameter sheets |
| `--proteinBinders` | `-pb` | `0` | Number of additional types of splicing factors to add to the simulation |

Example

```
python exportGeneSheet.py -g SRSF6 -gd /Data2/zane/human_genome/GCF_000001405.38/ -exp site -pb 3 -od /home/zane/Models/
```

### ```runGene.py```: Running a simulation for a single gene

| Argument | Short | Default | Description |
|-------------|------|---------|-------------|
| `--geneID` | `-g` | - | Name of gene from human genome |
| `--genomeDir` | `-gd` | - | Directory containing extracted human genome files (e.g. SRSF6.xlsx) |
| `--outputdir` | `-od` | `"./"` | Directory to save trajectory files |
| `--replicates` | `-r` | `1` | Number of replicate mRNA to simulate |
| `--coTrsc` | - | `True` | Use ```--coTrsc``` to include cotranscriptional splicing or ```--no-coTrsc``` to exclude |
| `--simTime` | `-t` | `60` | Amount of biological time to simulate per pre-mRNA in seconds. If coTrsc is True, this value will be added to the amount of time required to transcribe the gene. |
| `--writeInterval` | `-wi` | `1.0` | Frequency at which the simulation will record its state in biological seconds. |
| `--paramFile` | `-pf` | `None` | User input kinetic parameter file in the format exported by ```exportGeneSheet.py``` |
| `--paramSet` | `-ps` | `None` | If multiple similar parameter sets are being tested, the user can create individual parameter files named like ```parameters_PS.xlsx``` and this variable sets the value for ```PS``` |
| `--paramDir` | `-pd` | `None` | Directory where parameter file for ```--paramSet``` is located |
| `--outputFile` | `-of` | `None` | User can provide custom names for ```.lm``` trajectory files. Otherwise isosplicer uses its own naming scheme. |
| `--kinetics` | `-k` | `"global"` | Use ```global``` for unioform kinetics or ```site``` for site-specific kinetics |

Example

```
python runGene.py -r 10 -g SRSF6 -od ./ -t 600 -wi 60 -gd /Data2/zane/human_genome/GCF_000001405.38/ -pf ./SRSF6_parameters.xlsx -of ./SRSF6_test.lm -k site
```
