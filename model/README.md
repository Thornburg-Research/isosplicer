# isosplicer
Module to simulate mRNA splicing reaction dynamics for the human genome.

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
| `--outputdir` | `-od` | `"./"` | Whether the file contains a header row |
| `--export` | `-exp` | `"global"` | Export global (global) or global+site-specific (site) parameter sheets |
| `--proteinBinders` | `-pb` | `0` | Number of additional types of splicing factors to add to the simulation |

Example

```
python exportGeneSheet.py -g SRSF6 -gd /Data2/zane/human_genome/GCF_000001405.38/ -exp site -pb 3 -od /home/zane/Models/
```

### ```runGene.py```: Running a simulation for a single gene

| Argument | Short | Default | Description |
|----------|------|---------|-------------|
| `--geneID` | `-g` | - | Name of gene from human genome |
| `--genomeDir` | `-gd` | - | Directory containing extracted human genome files (e.g. SRSF6.xlsx) |
| `--outputdir` | `-od` | `"./"` | Whether the file contains a header row |
| `--outputdir` | `-od` | `"./"` | Whether the file contains a header row |
| `--export` | `-exp` | `"global"` | Export global (global) or global+site-specific (site) parameter sheets |
| `--proteinBinders` | `-pb` | `0` | Number of additional types of splicing factors to add to the simulation |

Example

```
python exportGeneSheet.py -g SRSF6 -gd /Data2/zane/human_genome/GCF_000001405.38/ -exp site -pb 3 -od /home/zane/Models/
```
