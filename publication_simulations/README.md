# Publication Simulation Executables

The following executables can be used with the input files provided here (along with the pre-processed human genome files) to run the simulations performed in the publication. 

Note that you might need to modify some of the paths, for example the path to your pre-processed human genome files (```-gd```) or the path to the gene lists (```-gl```).

## Multi-Isoform Genes Co-transcriptional Splicing

```
python runGenes.py -r 10 -th 20 -od ./Multi_Isoform_Genes_CoTrsc/ --coTrsc -t 600 -wi 60.0 -gl ./multiIsoformGenes.txt -gd /PATH/GCF_000001405.38/
```

## Multi-Isoform Genes Post-transcriptional Splicing

```
python runGenes.py -r 10 -th 20 -od ./Multi_Isoform_Genes_PostTrsc/ --no-coTrsc -t 600 -wi 60.0 -gl ./multiIsoformGenes.txt -gd /PATH/GCF_000001405.38/
```

## SRSF6 MCF7 Gradient Normoxia

```
python runGenes.py -r 1000 -th 3 -od ./SRSF6_MCF7_Grad_Norm/ --coTrsc -t 1200 -wi 60.0 -gl ./SRSF6.txt -gd /PATH/GCF_000001405.38/ -pd ./params_SRSF6_MCF7_Grad_Norm/ -ps 3 -k site
```

## SRSF6 MCF7 Gradient Hypoxia

```
python runGenes.py -r 1000 -th 3 -od ./SRSF6_MCF7_Grad_Hyp/ --coTrsc -t 1200 -wi 60.0 -gl ./SRSF6.txt -gd /PATH/GCF_000001405.38/ -pd ./params_SRSF6_MCF7_Grad_Hyp/ -ps 3 -k site
```

## SRSF6 MCF7 Phosphorylation Normoxia

```
python runGenes.py -r 1000 -th 3 -od ./SRSF6_MCF7_Phos_Norm/ --coTrsc -t 1200 -wi 60.0 -gl ./SRSF6.txt -gd /PATH/GCF_000001405.38/ -pd ./params_SRSF6_MCF7_Phos_Norm/ -ps 3 -k site
```

## SRSF6 MCF7 Phosphorylation Hypoxia

```
python runGenes.py -r 1000 -th 3 -od ./SRSF6_MCF7_Phos_Hyp/ --coTrsc -t 1200 -wi 60.0 -gl ./SRSF6.txt -gd /PATH/GCF_000001405.38/ -pd ./params_SRSF6_MCF7_Phos_Hyp/ -ps 3 -k site
```

## SRSF6 Multi Cell Normoxia

```
python runGenes.py -r 1000 -th 3 -od ./SRSF6_Multi_Cell_Norm/ --coTrsc -t 1200 -wi 60.0 -gl ./SRSF6.txt -gd /PATH/GCF_000001405.38/ -pd ./params_SRSF6_Multi_Cell_Norm/ -ps 3 -k site
```

## SRSF6 Multi Cell Hypoxia

```
python runGenes.py -r 1000 -th 3 -od ./SRSF6_Multi_Cell_Hyp/ --coTrsc -t 1200 -wi 60.0 -gl ./SRSF6.txt -gd /PATH/GCF_000001405.38/ -pd ./params_SRSF6_Multi_Cell_Hyp/ -ps 3 -k site
```
