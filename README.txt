The scripts in the repository 'mavrich_hatfull_nature_micro_2017' contain the analyses 
used in the publication Mavrich & Hatfull, Nature Microbiology, 2017. Below is a brief
description of each file in the repository. 

Notes: 
- some scripts may have tools not utilized for this particular publication.
- the format of SuppData2 file is a compilation of data from different parts of these
analyses, and does not itself serve as input data for the scripts.



1. parse_genbank_file_data.py

Python3 script to extract data from genbank-formatted files, such as CDS features, host data, and 
viral taxonomy data.



2. retrieve_ncbi_taxonomy.py

Python3 script to retrieve taxonomic data from NCBI, such as for viral hosts parsed from the viral 
genbank-formatted files.



3. compute_mash_distance.sh

Bash script to compute nucleotide distances between whole genomes using Mash.



4. process_mash_data.py

Python3 script to process the output from compute_mash_distance.sh, in which duplicated
data is removed and the reference and query file paths are trimmed.



5. analyze_pham_data.py

Python3 script to perform several computations based on Phamerator-derived gene family
(pham) data, such as pairwise gene content dissimilarity.



6. compute_gene_specific_mash.py

Python3 script to compute nucleotide distances between homologous (shared) and 
non-homologous (non-shared) nucleotide sequences using Mash.



7. analyze_mash_network.py

Python3 script to predict evolutionary mode for each phage using mash and pham data.



8. analyze_vog_data.py

Python3 script to compute pairwise gene content dissimilarity using VOG data.



9. edit_phylogenetic_tree.py

Python3 script to manipulate phylogenetic trees, such as naming internal nodes to be
compatible with Count output data.



10. analyze_count_data.py

Python3 script to analyze gene family gains and losses from Count based on a
phylogenetic tree.



11. analyze_all_data.R

RStudio script to compile all datasets, compute misc. analyses (some of which are
then outputted and used in downstream python scripts), and visualize/plot data.


