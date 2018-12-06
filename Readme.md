## CoMW

Comparative Metatranscriptomics Workflow is a standardized and validated workflow to functionally classify quality filtered mRNA reads from metatranscriptomic or total RNA studiesgenerated using HiSeq, MiSeq. CoMW is used for classification of these reads using assembled contigs to the reference databases provided and cited. 

If you use CoMW in your research, please cite:


For queries use issues or contact : mzanwar@envs.au.dk

# Installation

1. Download the latest distribution and expand it using unzip

	OR 

2. Download the development version:
   > git clone https://github.com/Zohaib509/CoMW.git



# Scripts

These scripts are written in Python and detailed parameters and dependdencies are along with usage examples are given in CoMW user manual CoMW_Manual.pdf

1. map_reads_to_contigs.py

2. filter_table_by_abundance.py [Optional]

3. align_contigs_to_database.py

4. parse_sword.py

5. map_orthologs_to_count_table.py

6. annotate_count_table.py



# Utils

Utils are small snippets written in R or bash to assist main scripts.
They must be present in utils folder in order to be accessible


# Databases 

Databases include fasta files and annotations that are also available in their respective developer's websites but are also collected here under same licence
Depending upon your usage, please cite the efforts of these databases that were developed by these groups.

1. The M5nr: a novel non-redundant database containing protein sequences and annotations from multiple sources and associated tools
Wilke, A., Harrison, T., Wilkening, J., Field, D., Glass, E.M., Kyrpides, N., Mavrommatis, K. and Meyer, F., 2012. The M5nr: a novel non-redundant database containing protein sequences and annotations from multiple sources and associated tools. BMC bioinformatics, 13(1), p.141.


2. Carbohydrate Active Enzymes database) and URL (http://www.cazy.org/) and cite :
Lombard, V., Golaconda Ramulu, H., Drula, E., Coutinho, P.M. and Henrissat, B., 2013. The carbohydrate-active enzymes database (CAZy) in 2013. Nucleic acids research, 42(D1), pp.D490-D495. 


3. NCycDB: a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes
Tu, Q., Lin, L., Cheng, L., Deng, Y., He, Z. and Wren, J., 2018. NCycDB: a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes. Bioinformatics, 1, p.9.


# Acknowledgment

This work was supported by a grant from the European Commission’s Marie Sklowdowska Curie Actions program MicroArctic-ITN under project number 675546.
