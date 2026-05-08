# Human Genome Preprocessing Notebooks
To run isosplicaer, sequence data for the human genome is needed. Here are the jupyter notebooks to process the human genome into the format that can be read by the model.

The scripts here were designed for the following human genome annotation: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.38/

Download the genome assembly, selecting the annotation features (gtf). This will download a zip file. Once you unzip the file, locate the following file: ```genomic.gtf```

To run, create a conda or python environment (python >=3.10).

Activate your environment and install the following dependencies:
```
pip install pyranges
pip install openpyxl
pip install jupyter
pip install numpy
pip install matplotlib
```

First, run the notebook ```Extract_Chromo_from_GTF.ipynb``` to isolate gene data into individual feature files. The ```genomic.gtf``` file is the input, change the input file path to your copy of the file. Also change the output path as desired.

Second, run ```Extract_Exons_from_GTF.ipynb``` to extract the individual gene exon data. The input is the output from ```Extract_Chromo_from_GTF.ipynb```, so change input path as needed. Also feel free to change output path as desired.

The output from ```Extract_Exons_from_GTF.ipynb``` are the files referenced for the splicing model.
