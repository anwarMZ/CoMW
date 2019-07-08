# CoMW

Comparative Metatranscriptomics Workflow is a standardized and validated workflow to functionally classify quality filtered mRNA reads from metatranscriptomic or total RNA studies generated using NGS short reads. CoMW is used for classification of these reads using assembled contigs to the reference databases provided and cited. 

If you use CoMW in your research, please cite:

Anwar MZ, Lanzen A, Bang-Andreasen T, and Jacobsen CS. To assemble or not to resemble – A validated Comparative Metatranscriptomics Workflow (CoMW), _GigaScience_,[Accepted, In Press]


For queries or issues please contact : mzanwar@envs.au.dk

## Installation

1. Download the latest distribution and expand it using unzip

	OR 

2. Download the development version:
   > git clone https://github.com/anwarMZ/CoMW.git

### Create an environment using anaconda, If you do not have anaconda installed, use [Anaconda installer link](https://docs.anaconda.com/anaconda/install/linux/)
> cd CoMW <br />
> conda env create -f ./CoMW.yml <br />
> source activate CoMW<br />

### Run install.sh file to download and make sword and databases to be used in CoMW
> bash ./test.sh

### Now run python scripts e.g.
> python scripts/assemble_reads.py -h 

## Scripts

These scripts are written in Python and detailed parameters and dependdencies are along with usage examples are given in CoMW user manual CoMW_Manual.pdf

1. assemble_reads.py

2. filter_ncRNAs.py

3. map_reads_to_contigs.py

4. filter_table_by_abundance.py [Optional]

5. align_contigs_to_database.py

6. parse_sword.py

7. map_orthologs_to_count_table.py

8. annotate_count_table.py


## Utils

Utils are small snippets written in R or bash to assist main scripts.
They must be present in utils folder in order to be accessible


## Databases 

Databases include fasta files and annotations that are also available in their respective developer's websites but are also collected here under same licence
Depending upon your usage, please cite the efforts of these databases that were developed by these groups.

1. The M5nr: a novel non-redundant database containing protein sequences and annotations from multiple sources and associated tools
Wilke, A., Harrison, T., Wilkening, J., Field, D., Glass, E.M., Kyrpides, N., Mavrommatis, K. and Meyer, F., 2012. The M5nr: a novel non-redundant database containing protein sequences and annotations from multiple sources and associated tools. BMC bioinformatics, 13(1), p.141.


2. Carbohydrate Active Enzymes database) and URL (http://www.cazy.org/) and cite :
Lombard, V., Golaconda Ramulu, H., Drula, E., Coutinho, P.M. and Henrissat, B., 2013. The carbohydrate-active enzymes database (CAZy) in 2013. Nucleic acids research, 42(D1), pp.D490-D495. 


3. NCycDB: a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes
Tu, Q., Lin, L., Cheng, L., Deng, Y., He, Z. and Wren, J., 2018. NCycDB: a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes. Bioinformatics, 1, p.9.


## Acknowledgment

This work was supported by a grant from the European Commission’s Marie Sklowdowska Curie Actions program MicroArctic-ITN under project number 675546.
