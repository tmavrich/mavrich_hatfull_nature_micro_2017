The repository 'mavrich_hatfull_nature_micro_2017' contain the scripts and input data
used for the analyses in the publication Mavrich & Hatfull, Nature Microbiology, 2017. 
Below is a brief description of each file in the repository. 


Scripts folder:

Notes: 
- some scripts may have tools not utilized for this particular publication.
- the format of SuppData2 file is a compilation of data from different parts of these
analyses, and does not itself serve as input data for the scripts.



1. parse_genbank_file_data.py

Description: Python3 script to extract data from genbank-formatted files, 
such as CDS features, host data, and viral taxonomy data.




2. retrieve_ncbi_taxonomy.py

Description: Python3 script to retrieve taxonomic data from NCBI, 
such as for viral hosts parsed from the viral genbank-formatted files.



3. compute_mash_distance.sh

Description: Bash script to compute nucleotide distances between whole genomes using Mash.



4. process_mash_data.py

Description: Python3 script to process the output from compute_mash_distance.sh, 
in which duplicated data is removed and the reference and query file paths are trimmed.



5. analyze_pham_data.py

Description: Python3 script to perform several computations based on Phamerator-derived 
gene family (pham) data, such as pairwise gene content dissimilarity.



6. compute_gene_specific_mash.py

Description: Python3 script to compute nucleotide distances between homologous (shared) 
and non-homologous (non-shared) nucleotide sequences using Mash. It does not compute all
pairwise genome comparisons, but only those specified in a required input file.



7. analyze_mash_network.py

Description: Python3 script to predict evolutionary mode for each phage using mash and
pham data.



8. analyze_vog_data.py

Description: Python3 script to compute pairwise gene content dissimilarity using VOG data.



9. edit_phylogenetic_tree.py

Description: Python3 script to manipulate phylogenetic tree data, such as naming 
internal nodes to be compatible with Count output data.



10. analyze_all_data.R

Description: RStudio script to compile all datasets, compute misc. analyses 
(some of which are then outputted and used in downstream python scripts), 
and visualize/plot data.











Input file folder:

Below is a list of the data files used in the R analysis and a brief description 
of how they were generated with the scripts listed above:



1. ani_79_data.csv

Description: Average nucleotide identity of 79 phage genomes used to optimize Mash
parameters. 

Data generation: Output from DNAMaster, manually edited in Excel.



2. ani_cluster_data.csv

Description: Average nucleotide identify of phages from specific clusters to assess
genomic similarities of HGCF and LGCF phages. 

Data generation: Output from DNAMaster, manually edited in Excel.



3. misc_clusters_gc.csv

Description: % GC data for several clusters.

Data generation: SQL query of Phamerator database, manually edited in Excel.



4. misc_clusters_genometrics.csv

Description: Misc. genometric data for several clusters.

Data generation: SQL query of Phamerator database, manually edited in Excel.



5. phage_host_data.csv

Description: Phage and host metadata.

Data generation: Retrieved from misc. sources, manually edited in Excel.



6. phylogeny_data.csv

Description: Pairwise branch length distances computed from phylogenetic tree.

Data generation: Seaview output, manually edited in Excel.



7. processed_mash_output.csv

Description: Mash distance data.

Data generation: 

Run process_mash_data.py script on:
    a.  output file from compute_mash_distance.sh script. To obtain mash output file, 
        run compute_mash_distance.sh on folder of genome fasta files.



8. pairwise_pham_proportions.csv

Description: Pairwise gene content dissimilarity data using phams created 
from Phamerator.

Data generation:

Run analyze_pham_data.py script on:
    a. genome file generated from SQL query of Phamerator database, manually edited in Excel.
    b. gene file generated from SQL query of Phamerator database, manually edited in Excel.



9. pairwise_vog_proportions.csv

Description: Pairwise gene content dissimilarity data using VOGs created
from the online pVOG database.

Data generation: 

Run analyze_vog_data.py script on:
    a. list of accessions of phages that need to be analyzed
    b. VOG data file downloaded from online database



10. actino_only_pairwise_pham_proportions.csv

Description: Pairwise gene content dissimilarity data using phams created 
from Phamerator for only Actinobacteriophages and includes pham function data.

Data generation: 

Run analyze_pham_data.py script on:
    a. genome file generated from SQL query of Phamerator database, manually edited in Excel.
    b. gene file generated from SQL query of Phamerator database, manually edited in Excel.



11. mode_prediction.csv

Description: Predicted evolutionary mode (HGCF, LGCF) of all phages in the analysis.

Data generation: 

Run analyze_mash_network.py script on:
    a. list of phage names of interest, manually edited in Excel.
    b. mash distance and gene content dissimilarity data generated from analyze_all_data.R script.



12. gene_specific_mash_data.csv

Description: Pairwise mash distances of shared (homologous) and unshared (non-homologous)
genes.

Data generation: 

Run compute_gene_specific_mash.py script on:
    a. folder of genome fasta files.
    b. gene file generated from SQL query of Phamerator database, manually edited in Excel.
    c. list of specific pairwise comparisons to analyze generated from analyze_all_data.R
        script.


