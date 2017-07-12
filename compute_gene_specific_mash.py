#Python3 script to compute mash distances based on gene-specific sequences
#Travis Mavrich
#20170105



#Import requisite modules
import time, sys, os, csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import subprocess

#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    fasta_file_dir = sys.argv[1]
    gene_file = sys.argv[2]
    comparison_file = sys.argv[3]

except:


    print("\n\n\
        This is a Python3 script to compute mash distances based on gene-specific sequences\n\
        Execute script in any working directory\n\
        Script requires three arguments:\n\n\
        First argument: fasta file directory file (fasta-formatted)\n\
            \n\
        Second argument: gene file (csv-formatted)\n\
            0 PhageID\n\
            1 GeneID\n\
            2 Gene Start (1-based indexing)\n\
            3 Gene Stop (1-based indexing)\n\
            4 Pham\n\
            \n\
        Third argument: comparison list file (csv-formatted)\n\
            0 PhageID 1\n\
            1 PhageID 2\n\
            \n\
            \n")


        
    sys.exit(1)



#Expand home and working directory
home_dir = os.path.expanduser('~')
working_dir = os.path.abspath('.')






#Verify the fasta file folder exists


#Expand the path if it references the home directory
if fasta_file_dir[0] == "~":
    fasta_file_dir = home_dir + fasta_file_dir[1:]

if fasta_file_dir[0] != "/":
    fasta_file_dir = working_dir + '/' + fasta_file_dir



#Expand the path, to make sure it is a complete directory path (in case user inputted path with './path/to/folder')
fasta_file_dir = os.path.abspath(fasta_file_dir)

if fasta_file_dir[-1] != "/":
    fasta_file_dir = fasta_file_dir + '/'


if os.path.isdir(fasta_file_dir) == False:
    print("\n\nInvalid input for fasta folder.\n\n")
    sys.exit(1)








#Verify the gene data file exists

#Expand the path if it references the home directory
if gene_file[0] == "~":
    gene_file = home_dir + gene_file[1:]

if gene_file[0] != "/":
    gene_file = working_dir + '/' + gene_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
gene_file = os.path.abspath(gene_file)

if os.path.exists(gene_file) == False:
    print("\n\nInvalid gene data file path.\n\n")
    sys.exit(1)



#Verify the comparison data file exists

#Expand the path if it references the home directory
if comparison_file[0] == "~":
    comparison_file = home_dir + comparison_file[1:]

if comparison_file[0] != "/":
    comparison_file = working_dir + '/' + comparison_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
comparison_file = os.path.abspath(comparison_file)

if os.path.exists(comparison_file) == False:
    print("\n\nInvalid comparison list file path.\n\n")
    sys.exit(1)





def make_fasta(seq_id,sequence,fasta_file_name):
    output_handle = open(fasta_file_name,'w')
    output_handle.write(">%s" % seq_id)
    output_handle.write('\n' + sequence)
    output_handle.close()




def compute_pairwise_mash(fasta1_filename,fasta2_filename):

    sketch_size = 25000
    kmer_size = 15

    #Compute sketches of each shared sequence
    command_string = "mash sketch -s %s -k %s %s" % (sketch_size,kmer_size,fasta1_filename)
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list)

    command_string = "mash sketch -s %s -k %s %s" % (sketch_size,kmer_size,fasta2_filename)
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list)

    #Compute mash distance of shared sequences
    fasta1_sketch_file_name = fasta1_filename + '.msh'
    fasta2_sketch_file_name = fasta2_filename + '.msh'

    mash_output_file = 'mash_output.txt'
    mash_output_file_handle = open(mash_output_file,'w')
    command_string = "mash dist %s %s" % (fasta1_sketch_file_name,fasta2_sketch_file_name)
    command_list = command_string.split(" ")
    proc = subprocess.check_call(command_list,stdout = mash_output_file_handle)
    mash_output_file_handle.close()


    #Import mash distance results
    mash_input_handle = open(mash_output_file,'r')
    mash_input_reader = csv.reader(mash_input_handle,delimiter = '\t')

    #Mash output structure:
    #0 reference
    #1 query
    #2 mash distance
    #3 mash p-value
    #4 kmer count
    mash_input_list = []


    for row in mash_input_reader:
        row[2] = float(row[2])
        row[3] = float(row[3])
        mash_input_list.append(row)
    mash_input_handle.close()

    #Remove all files
    #input('About to delete mash files')
    os.remove(fasta1_filename)
    os.remove(fasta2_filename)
    os.remove(fasta1_sketch_file_name)
    os.remove(fasta2_sketch_file_name)
    os.remove(mash_output_file)
    

    return mash_input_list






date = time.strftime("%Y%m%d")




#Import raw data
print("Importing data")
comparison_data_description = input("Provide description of the comparison data set: ")




#Import each fasta file, and save genome sequences in a dictionary
#The script only expects one genome sequence per file
fasta_dict = {}
files =  [X for X in os.listdir(fasta_file_dir) if os.path.isfile(os.path.join(fasta_file_dir,X))]

for filename in files:

    #The filename should contain the phageID that will be used to match to other data
    fasta_phageID = filename.split('.')[0]
    seqFile = fasta_file_dir + filename
    
    if fasta_phageID in fasta_dict:
        print("Error: %s is a duplicate fasta phageID. Skipping %s fasta record." % (fasta_phageID,filename))
        input("Press ENTER to continue")
        continue
    
    #Create a fasta record and store in the fasta dictionary
    fasta_record = SeqIO.read(seqFile, "fasta")
    fasta_dict[fasta_phageID] = fasta_record




#Import gene data file
file_object = open(gene_file,"r")
file_reader = csv.reader(file_object)
gene_data_list = []
for row in file_reader:
    gene_data_list.append(row)
file_object.close()




#Import list of comparisons
file_object = open(comparison_file,"r")
file_reader = csv.reader(file_object)
comparison_data_list = []
for row in file_reader:
    comparison_data_list.append(row)
file_object.close()





#Create output directory to store all analyses
new_dir = date + "_gene_specific_mash_analysis"
os.mkdir(new_dir)
os.chdir(new_dir)






#Create dictionaries, sets, and lists to perform computations


#Create phageID set from comparison list
#This will serve as a record of all phages that will be used in the analysis.
#This will be used to:
#1. verify that all comparison phageIDs match to fasta record phageIDs and gene data phageIDs
#2. determine which genomes from the fasta records actually will be used to extract gene sequences (in case there are fasta records that are unused in the analysis)
print("Creating comparison phageID set")
comparison_phageID_set = set()
for comparison_data in comparison_data_list:
    comparison_phageID_set.add(comparison_data[0])
    comparison_phageID_set.add(comparison_data[1])



#Create phage-gene dictionary and phage-pham dictionary
print("Creating phage-gene and phage-pham dictionaries")

#First create the set of all phageIDs in the gene data list
gene_data_phageID_set = set()
for gene_data in gene_data_list:
    gene_data_phageID_set.add(gene_data[0])

#Now, using all unique gene data phageIDs, create a gene data dictionary
#Key = phageID
#Value = list of gene data
#0 = phageID
#1 = GeneID
#2 = Gene start
#3 = Gene stop
#4 = Pham


#And create phage-pham dictionary
#Key = phageID
#Value = set of phams


phage_gene_dict = {}
phage_pham_dict = {}
for gene_data_phageID in gene_data_phageID_set:
    matching_gene_list = []
    matching_pham_set = set()

    for gene_data in gene_data_list:
        if gene_data[0] == gene_data_phageID:
            matching_gene_list.append(gene_data)
            matching_pham_set.add(gene_data[4])

    phage_gene_dict[gene_data_phageID] = matching_gene_list
    phage_pham_dict[gene_data_phageID] = matching_pham_set










#Now that all phageID sets have been created, verify that there is gene data and genome sequence data for each comparison

diff_set1 = comparison_phageID_set - gene_data_phageID_set
if len(diff_set1) > 0:
    print("Warning: there are %s phageIDs in the comparison data that have no gene data." % len(diff_set1))
    print(diff_set1)
    input("Press ENTER to continue")


diff_set2 = comparison_phageID_set - fasta_dict.keys()
if len(diff_set2) > 0:
    print("Warning: there are %s phageIDs in the comparison data that have no genome sequence." % len(diff_set2))
    print(diff_set2)
    input("Press ENTER to continue")









#Create phage-gene sequence dictionary for only the genomes needed in the comparison list
#Key = phageID
#Value = list of gene sequence data
#0 = phageID
#1 = GeneID
#2 = Pham
#3 = Gene sequence
print("Creating phage-gene sequence dictionary")
phage_gene_sequence_dict = {}

for comparison_phageID in comparison_phageID_set:
    matching_gene_sequence_list = []
    fasta_seq = fasta_dict[comparison_phageID]

    for gene_data in phage_gene_dict[comparison_phageID]:

        #Extract gene sequence from fasta record
        #Slicing includes the first coordinate and excludes the last
        #The coordinates in the input file are expected to be 1-based indexing of the first and last nucleotides
        #Therefore, the start coordinate needs to be reduced by 1, and the slice should extract the expected sequence
        gene_sequence = fasta_seq.seq[int(gene_data[2]) - 1:int(gene_data[3])]        
        gene_sequence = str(gene_sequence).upper() #Convert to string to get rid of seq object
        
        #Print the details of the sequence extraction to ensure it's working correctly
        #print(gene_data[0])
        #print(gene_data[1])
        #print(int(gene_data[2]))
        #print(int(gene_data[3]))
        #print(len(gene_sequence))
        #print(gene_sequence)
        #input()
        
        #Now create a new list of gene data, including the sequence and append to list of all gene data for that phage
        gene_sequence_data = []        
        gene_sequence_data.append(gene_data[0])
        gene_sequence_data.append(gene_data[1])
        gene_sequence_data.append(gene_data[4])
        gene_sequence_data.append(gene_sequence)
        matching_gene_sequence_list.append(gene_sequence_data)

    phage_gene_sequence_dict[comparison_phageID] = matching_gene_sequence_list













#Now, for each comparison to analyze:
#1: create sets of shared and unshared phams
#2: create concatenated nucleotide sequences for all, shared, and unshared gene sequences for each phage in the comparison
#3: compute mash distances of all, shared, and unshared concatenated sequences
print("Computing all, shared, and unshared mash distances")

#This list will stored results for shared data
#0 phage1_phage2
#1 phage1
#2 phage2
#3 phage1_phage2 # shared phams
#4 phage1_phage2 gene content dissimilarity (general index)
#5 phage1_phage2 gene content dissimilarity (jaccard index)
#6 phage1_phage2 all genes mash distance
#7 phage1_phage2 all genes mash p-value
#8 phage1_phage2 all genes mash kmer count
#9 phage1_phage2 shared genes mash distance
#10 phage1_phage2 shared genes mash p-value
#11 phage1_phage2 shared genes mash kmer count
#12 phage1_phage2 unshared genes mash distance
#13 phage1_phage2 unshared genes mash p-value
#14 phage1_phage2 unshared genes mash kmer count
#15 phage1_phage2 shared-unshared mash distance
#16 phage1_phage2 shared-unshared mash p-value
#17 phage1_phage2 shared-unshared mash kmer count
#18 phage1_phage2 total length of combined shared sequence
#19 phage1_phage2 total length of combined unshared sequence
#20 phage1_phage2 combined shared gene GC content
#21 phage1_phage2 combined unshared gene GC content

#22 phage1 # unshared phams
#23 phage1 # all phams
#24 phage1 # shared genes
#25 phage1 # unshared genes
#26 phage1 # all genes
#27 phage1 average length of all genes
#28 phage1 average length of shared genes
#29 phage1 average length of unshared genes
#30 phage1 total length of all genes
#31 phage1 total length of shared genes
#32 phage1 total length of unshared genes
#33 phage1 all genes GC content
#34 phage1 shared genes GC content
#35 phage1 unshared genes GC content
    
#36 phage2 # unshared phams
#37 phage2 # all phams
#38 phage2 # shared genes
#39 phage2 # unshared genes
#40 phage2 # all genes
#41 phage2 average length of all genes
#42 phage2 average length of shared genes
#43 phage2 average length of unshared genes
#44 phage2 total length of all genes
#45 phage2 total length of shared genes
#46 phage2 total length of unshared genes
#47 phage2 all genes GC content
#48 phage2 shared genes GC content
#49 phage2 unshared genes GC content



shared_unshared_data_header = []
shared_unshared_data_header.append(['Shared and unshared pham and gene data'])
shared_unshared_data_header.append([comparison_data_description])
shared_unshared_data_header.append([date])

shared_unshared_data_columns = ['phage1_phage2',\
                                'phage1',\
                                'phage2',\
                                'phage1_phage2 # shared phams',\
                                'phage1_phage2 gene content dissimilarity (general index)',\
                                'phage1_phage2 gene content dissimilarity (jaccard index)',\
                                'phage1_phage2 all genes mash distance',\
                                'phage1_phage2 all genes mash p-value',\
                                'phage1_phage2 all genes mash kmer count',\
                                'phage1_phage2 shared genes mash distance',\
                                'phage1_phage2 shared genes mash p-value',\
                                'phage1_phage2 shared genes mash kmer count',\
                                'phage1_phage2 unshared genes mash distance',\
                                'phage1_phage2 unshared genes mash p-value',\
                                'phage1_phage2 unshared genes mash kmer count',\
                                'phage1_phage2 shared-unshared mash distance',\
                                'phage1_phage2 shared-unshared mash p-value',\
                                'phage1_phage2 shared-unshared mash kmer count',\
                                'phage1_phage2 total length of combined shared sequence',\
                                'phage1_phage2 total length of combined unshared sequence',\
                                'phage1_phage2 combined shared gene GC content',\
                                'phage1_phage2 combined unshared gene GC content',\
                                'phage1 # unshared phams',\
                                'phage1 # all phams',\
                                'phage1 # all genes',\
                                'phage1 # shared genes',\
                                'phage1 # unshared genes',\
                                'phage1 average length of all genes',\
                                'phage1 average length of shared genes',\
                                'phage1 average length of unshared genes',\
                                'phage1 total length of all genes',\
                                'phage1 total length of shared genes',\
                                'phage1 total length of unshared genes',\
                                'phage1 all genes GC content',\
                                'phage1 shared genes GC content',\
                                'phage1 unshared genes GC content',\
                                'phage2 # unshared phams',\
                                'phage2 # all phams',\
                                'phage2 # all genes',\
                                'phage2 # shared genes',\
                                'phage2 # unshared genes',\
                                'phage2 average length of all genes',\
                                'phage2 average length of shared genes',\
                                'phage2 average length of unshared genes',\
                                'phage2 total length of all genes',\
                                'phage2 total length of shared genes',\
                                'phage2 total length of unshared genes',\
                                'phage2 all genes GC content',\
                                'phage2 shared genes GC content',\
                                'phage2 unshared genes GC content']


output_file_name =  new_dir + '.csv'
output_csvfile = open(output_file_name,'w')
output_writer = csv.writer(output_csvfile)
for row in shared_unshared_data_header:
    output_writer.writerow(row)
output_writer.writerow(shared_unshared_data_columns)


temp_dir = "temp_dir"
os.mkdir(temp_dir)
os.chdir(temp_dir)


comparison_count = 0
for comparison_data in comparison_data_list:

    shared_unshared_gene_data_results = []

    phage1 = comparison_data[0]
    phage2 = comparison_data[1]
    phage1_phage2 = phage1 + "_" + phage2

    comparison_count += 1
    print("Computing shared and unshared mash distances for comparison #%s" % comparison_count)
    #input("Press ENTER to continue")

    phage1_pham_set = phage_pham_dict[phage1]
    phage1_gene_sequence_list = phage_gene_sequence_dict[phage1]
    
    phage2_pham_set = phage_pham_dict[phage2]
    phage2_gene_sequence_list = phage_gene_sequence_dict[phage2]



    #If there are no phams in either genome, then skip this comparison
    if len(phage1_pham_set) == 0:
        print("Phage %s contains no phams." % phage1)
        continue

    if len(phage2_pham_set) == 0:
        print("Phage %s contains no phams." % phage2)
        continue
    
    #If both phages to compare are the same, then skip this genome
    if phage1 == phage2:
        print("Phages %s and %s are the same." % (phage1,phage2))
        continue



    phage1_unshared_pham_set = phage1_pham_set - phage2_pham_set
    phage2_unshared_pham_set = phage2_pham_set - phage1_pham_set
    
    intersection_set = phage1_pham_set & phage2_pham_set
    union_set = phage1_pham_set | phage2_pham_set

    phage1_unshared_size = len(phage1_unshared_pham_set)
    phage2_unshared_size = len(phage2_unshared_pham_set)
    phage1_phage2_shared_size = len(intersection_set)
    phage1_shared_proportion = len(intersection_set)/len(phage1_pham_set)
    phage2_shared_proportion = len(intersection_set)/len(phage2_pham_set)
    ave_proportion = (phage1_shared_proportion + phage2_shared_proportion)/2
    jaccard_similarity = len(intersection_set)/len(union_set)
    gene_content_dissimilarity_general = 1 - ave_proportion
    gene_content_dissimilarity_jaccard = 1 - jaccard_similarity


    if len(phage1_unshared_pham_set) + len(phage2_unshared_pham_set) + len(intersection_set) != len(union_set):
        print("Error in calculating shared and unshared pham sets for %s" % phage1_phage2)
        continue
    


    #Retrieve all genes for each group of shared and unshared phams, and concatenate into a single long sequence
    phage1_all_gene_sequence = ""
    phage2_all_gene_sequence = ""

    phage1_shared_gene_sequence = ""
    phage2_shared_gene_sequence = ""

    phage1_unshared_gene_sequence = ""    
    phage2_unshared_gene_sequence = ""

    phage1_all_gene_lengths = []
    phage2_all_gene_lengths = []

    phage1_shared_gene_lengths = []
    phage2_shared_gene_lengths = []

    phage1_unshared_gene_lengths = []
    phage2_unshared_gene_lengths = []

    phage1_all_gene_count = 0
    phage2_all_gene_count = 0

    phage1_shared_gene_count = 0
    phage2_shared_gene_count = 0

    phage1_unshared_gene_count = 0
    phage2_unshared_gene_count = 0




    for phage_gene_data in phage1_gene_sequence_list:

        phage1_all_gene_count += 1        
        phage1_all_gene_sequence = phage1_all_gene_sequence + phage_gene_data[3]
        phage1_all_gene_lengths.append(len(phage_gene_data[3]))

        
        if phage_gene_data[2] in intersection_set:
            phage1_shared_gene_count += 1
            phage1_shared_gene_lengths.append(len(phage_gene_data[3]))
            phage1_shared_gene_sequence = phage1_shared_gene_sequence + phage_gene_data[3]
    
        if phage_gene_data[2] in phage1_unshared_pham_set:
            phage1_unshared_gene_count += 1
            phage1_unshared_gene_lengths.append(len(phage_gene_data[3]))
            phage1_unshared_gene_sequence = phage1_unshared_gene_sequence + phage_gene_data[3]


    for phage_gene_data in phage2_gene_sequence_list:

        phage2_all_gene_count += 1        
        phage2_all_gene_sequence = phage2_all_gene_sequence + phage_gene_data[3]
        phage2_all_gene_lengths.append(len(phage_gene_data[3]))

        
        if phage_gene_data[2] in intersection_set:
            phage2_shared_gene_count += 1
            phage2_shared_gene_lengths.append(len(phage_gene_data[3]))
            phage2_shared_gene_sequence = phage2_shared_gene_sequence + phage_gene_data[3]
    
        if phage_gene_data[2] in phage2_unshared_pham_set:
            phage2_unshared_gene_count += 1
            phage2_unshared_gene_lengths.append(len(phage_gene_data[3]))
            phage2_unshared_gene_sequence = phage2_unshared_gene_sequence + phage_gene_data[3]



    phage1_all_genes_total_length = sum(phage1_all_gene_lengths)
    phage2_all_genes_total_length = sum(phage2_all_gene_lengths)

    if phage1_all_gene_count > 0:
        phage1_all_genes_ave_length = phage1_all_genes_total_length/len(phage1_all_gene_lengths)
        
        
        #Print details to ensure it is working correctly
        #print(phage1)
        #print(phage1_all_gene_lengths)
        #print(len(phage1_all_gene_lengths))
        #print(phage1_all_genes_total_length)
        #print(phage1_all_genes_ave_length)
        #input()
        
    else:
        phage1_all_genes_ave_length = 0

    if phage2_all_gene_count > 0:
        phage2_all_genes_ave_length = phage2_all_genes_total_length/len(phage2_all_gene_lengths)

        #Print details to ensure it is working correctly
        #print(phage2)
        #print(phage2_all_gene_lengths)
        #print(len(phage2_all_gene_lengths))
        #print(phage2_all_genes_total_length)
        #print(phage2_all_genes_ave_length)
        #input()


    else:
        phage2_all_genes_ave_length = 0



    phage1_shared_genes_total_length = sum(phage1_shared_gene_lengths)
    phage2_shared_genes_total_length = sum(phage2_shared_gene_lengths)

    if phage1_shared_gene_count > 0:
        phage1_shared_genes_ave_length = phage1_shared_genes_total_length/len(phage1_shared_gene_lengths)

        #Print details to ensure it is working correctly
        #print(phage1)
        #print(phage1_shared_gene_lengths)
        #print(len(phage1_shared_gene_lengths))
        #print(phage1_shared_genes_total_length)
        #print(phage1_shared_genes_ave_length)
        #input()


    else:
        phage1_shared_genes_ave_length = 0

    if phage2_shared_gene_count > 0:
        phage2_shared_genes_ave_length = phage2_shared_genes_total_length/len(phage2_shared_gene_lengths)

        #Print details to ensure it is working correctly
        #print(phage2)
        #print(phage2_shared_gene_lengths)
        #print(len(phage2_shared_gene_lengths))
        #print(phage2_shared_genes_total_length)
        #print(phage2_shared_genes_ave_length)
        #input()


    else:
        phage2_shared_genes_ave_length = 0






    phage1_unshared_genes_total_length = sum(phage1_unshared_gene_lengths)
    phage2_unshared_genes_total_length = sum(phage2_unshared_gene_lengths)


    if phage1_unshared_gene_count > 0:
        phage1_unshared_genes_ave_length = phage1_unshared_genes_total_length/len(phage1_unshared_gene_lengths)

        #Print details to ensure it is working correctly
        #print(phage1)
        #print(phage1_unshared_gene_lengths)
        #print(len(phage1_unshared_gene_lengths))
        #print(phage1_unshared_genes_total_length)
        #print(phage1_unshared_genes_ave_length)
        #input()



    else:
        phage1_unshared_genes_ave_length = 0

    if phage2_unshared_gene_count > 0:
        phage2_unshared_genes_ave_length = phage2_unshared_genes_total_length/len(phage2_unshared_gene_lengths)

        #Print details to ensure it is working correctly
        #print(phage2)
        #print(phage2_unshared_gene_lengths)
        #print(len(phage2_unshared_gene_lengths))
        #print(phage2_unshared_genes_total_length)
        #print(phage2_unshared_genes_ave_length)
        #input()


    else:
        phage2_unshared_genes_ave_length = 0



    phage1_phage2_shared_gene_sequence = phage1_shared_gene_sequence + phage2_shared_gene_sequence #combined sequence of shared genes from both genomes
    phage1_phage2_unshared_gene_sequence = phage1_unshared_gene_sequence + phage2_unshared_gene_sequence #combined sequence of unshared genes from both genomes
    phage1_phage2_shared_gene_sequence_total_length = len(phage1_phage2_shared_gene_sequence)
    phage1_phage2_unshared_gene_sequence_total_length = len(phage1_phage2_unshared_gene_sequence)


    #Compute GC content for all concatenated sequences
    #Be sure to test the length of the nucleotide sequence before computing GC%, since some shared or unshared sequences make have length of zero

    if len(phage1_all_gene_sequence) > 0:
        phage1_all_gene_GC_content = (phage1_all_gene_sequence.count('G') + phage1_all_gene_sequence.count('C'))/len(phage1_all_gene_sequence)
    else:
        phage1_all_gene_GC_content = 0

    if len(phage2_all_gene_sequence) > 0:
        phage2_all_gene_GC_content = (phage2_all_gene_sequence.count('G') + phage2_all_gene_sequence.count('C'))/len(phage2_all_gene_sequence)
    else:
        phage2_all_gene_GC_content = 0



    if len(phage1_shared_gene_sequence) > 0:
        phage1_shared_gene_GC_content = (phage1_shared_gene_sequence.count('G') + phage1_shared_gene_sequence.count('C'))/len(phage1_shared_gene_sequence)
    else:
        phage1_shared_gene_GC_content = 0

    if len(phage2_shared_gene_sequence) > 0:
        phage2_shared_gene_GC_content = (phage2_shared_gene_sequence.count('G') + phage2_shared_gene_sequence.count('C'))/len(phage2_shared_gene_sequence)
    else:
        phage2_shared_gene_GC_content = 0
    
    
    
    if len(phage1_unshared_gene_sequence) > 0:
        phage1_unshared_gene_GC_content = (phage1_unshared_gene_sequence.count('G') + phage1_unshared_gene_sequence.count('C'))/len(phage1_unshared_gene_sequence)
    else:
        phage1_unshared_gene_GC_content = 0
        
    if len(phage2_unshared_gene_sequence) > 0:
        phage2_unshared_gene_GC_content = (phage2_unshared_gene_sequence.count('G') + phage2_unshared_gene_sequence.count('C'))/len(phage2_unshared_gene_sequence)
    else:
        phage2_unshared_gene_GC_content = 0



    if len(phage1_phage2_shared_gene_sequence) > 0:
        phage1_phage2_shared_gene_GC_content = (phage1_phage2_shared_gene_sequence.count('G') + phage1_phage2_shared_gene_sequence.count('C'))/len(phage1_phage2_shared_gene_sequence)
    else:
        phage1_phage2_shared_gene_GC_content = 0

    if len(phage1_phage2_unshared_gene_sequence) > 0:
        phage1_phage2_unshared_gene_GC_content = (phage1_phage2_unshared_gene_sequence.count('G') + phage1_phage2_unshared_gene_sequence.count('C'))/len(phage1_phage2_unshared_gene_sequence)
    else:
        phage1_phage2_unshared_gene_GC_content = 0





    #Compute mash distances of all, shared, unshared, and combined shared-combined unshared gene sequences


    #Distance between whole genome sequences
    if (len(phage1_all_gene_sequence) > 0 and len(phage2_all_gene_sequence) > 0):

        #Create all sequence fasta files
        phage1_all_fasta_file_name =  'phage1_all_sequence.fasta'
        make_fasta(phage1,phage1_all_gene_sequence,phage1_all_fasta_file_name)

        phage2_all_fasta_file_name =  'phage2_all_sequence.fasta'
        make_fasta(phage2,phage2_all_gene_sequence,phage2_all_fasta_file_name)

        #Compute all sequence distance
        all_sequence_mash_result = compute_pairwise_mash(phage1_all_fasta_file_name,phage2_all_fasta_file_name)

    else:
        all_sequence_mash_result = [[phage1,phage2,float(1),float(1),'0/0']]



    #Distance between shared sequences
    if (len(phage1_shared_gene_sequence) > 0 and len(phage2_shared_gene_sequence) > 0):

        #Create shared sequence fasta files
        phage1_shared_fasta_file_name =  'phage1_shared_sequence.fasta'
        make_fasta(phage1,phage1_shared_gene_sequence,phage1_shared_fasta_file_name)

        phage2_shared_fasta_file_name =  'phage2_shared_sequence.fasta'
        make_fasta(phage2,phage2_shared_gene_sequence,phage2_shared_fasta_file_name)

        #Compute shared sequence distance
        shared_sequence_mash_result = compute_pairwise_mash(phage1_shared_fasta_file_name,phage2_shared_fasta_file_name)

    else:
        shared_sequence_mash_result = [[phage1,phage2,float(1),float(1),'0/0']]



    #Distance between unshared sequences
    if (len(phage1_unshared_gene_sequence) > 0 and len(phage2_unshared_gene_sequence) > 0):

        #Create unshared sequence fasta files
        phage1_unshared_fasta_file_name =  'phage1_unshared_sequence.fasta'
        make_fasta(phage1,phage1_unshared_gene_sequence,phage1_unshared_fasta_file_name)

        phage2_unshared_fasta_file_name =  'phage2_unshared_sequence.fasta'
        make_fasta(phage2,phage2_unshared_gene_sequence,phage2_unshared_fasta_file_name)

        #Compute unshared sequence distance
        unshared_sequence_mash_result = compute_pairwise_mash(phage1_unshared_fasta_file_name,phage2_unshared_fasta_file_name)

    else:
        unshared_sequence_mash_result = [[phage1,phage2,float(1),float(1),'0/0']]



    #Distance between combined sequences
    if (len(phage1_phage2_shared_gene_sequence) > 0 and len(phage1_phage2_unshared_gene_sequence) > 0):

        #Create combined shared and combined unshared sequence fasta files
        phage1_phage2_shared_fasta_file_name =  'phage1_phage2_shared_sequence.fasta'
        make_fasta(phage1_phage2,phage1_phage2_shared_gene_sequence,phage1_phage2_shared_fasta_file_name)

        phage1_phage2_unshared_fasta_file_name =  'phage1_phage2_unshared_sequence.fasta'
        make_fasta(phage1_phage2,phage1_phage2_unshared_gene_sequence,phage1_phage2_unshared_fasta_file_name)

        #Compute shared-unshared sequence distance
        shared_unshared_sequence_mash_result = compute_pairwise_mash(phage1_phage2_shared_fasta_file_name,phage1_phage2_unshared_fasta_file_name)

    else:
        shared_unshared_sequence_mash_result = [[phage1,phage2,float(1),float(1),'0/0']]




    #Prepare data for export
    shared_unshared_gene_data_results.append(phage1_phage2)
    shared_unshared_gene_data_results.append(phage1)
    shared_unshared_gene_data_results.append(phage2)    
   
    shared_unshared_gene_data_results.append(len(intersection_set))
    shared_unshared_gene_data_results.append(round(gene_content_dissimilarity_general,4))    
    shared_unshared_gene_data_results.append(round(gene_content_dissimilarity_jaccard,4))


    shared_unshared_gene_data_results.append(round(all_sequence_mash_result[0][2],4))
    shared_unshared_gene_data_results.append(round(all_sequence_mash_result[0][3],4))    
    shared_unshared_gene_data_results.append(all_sequence_mash_result[0][4])
    shared_unshared_gene_data_results.append(round(shared_sequence_mash_result[0][2],4))
    shared_unshared_gene_data_results.append(round(shared_sequence_mash_result[0][3],4))    
    shared_unshared_gene_data_results.append(shared_sequence_mash_result[0][4])
    shared_unshared_gene_data_results.append(round(unshared_sequence_mash_result[0][2],4))
    shared_unshared_gene_data_results.append(round(unshared_sequence_mash_result[0][3],4))    
    shared_unshared_gene_data_results.append(unshared_sequence_mash_result[0][4])

    shared_unshared_gene_data_results.append(round(shared_unshared_sequence_mash_result[0][2],4))
    shared_unshared_gene_data_results.append(round(shared_unshared_sequence_mash_result[0][3],4))    
    shared_unshared_gene_data_results.append(shared_unshared_sequence_mash_result[0][4])
    shared_unshared_gene_data_results.append(phage1_phage2_shared_gene_sequence_total_length)
    shared_unshared_gene_data_results.append(phage1_phage2_unshared_gene_sequence_total_length)
    shared_unshared_gene_data_results.append(phage1_phage2_shared_gene_GC_content)
    shared_unshared_gene_data_results.append(phage1_phage2_unshared_gene_GC_content)
    
    shared_unshared_gene_data_results.append(phage1_unshared_size)
    shared_unshared_gene_data_results.append(len(phage1_pham_set))
    shared_unshared_gene_data_results.append(phage1_all_gene_count)
    shared_unshared_gene_data_results.append(phage1_shared_gene_count)    
    shared_unshared_gene_data_results.append(phage1_unshared_gene_count)
    shared_unshared_gene_data_results.append(round(phage1_all_genes_ave_length,4))
    shared_unshared_gene_data_results.append(round(phage1_shared_genes_ave_length,4))
    shared_unshared_gene_data_results.append(round(phage1_unshared_genes_ave_length,4))    
    shared_unshared_gene_data_results.append(phage1_all_genes_total_length)
    shared_unshared_gene_data_results.append(phage1_shared_genes_total_length)
    shared_unshared_gene_data_results.append(phage1_unshared_genes_total_length)
    shared_unshared_gene_data_results.append(round(phage1_all_gene_GC_content,4))
    shared_unshared_gene_data_results.append(round(phage1_shared_gene_GC_content,4))
    shared_unshared_gene_data_results.append(round(phage1_unshared_gene_GC_content,4))

    shared_unshared_gene_data_results.append(phage2_unshared_size)
    shared_unshared_gene_data_results.append(len(phage2_pham_set))
    shared_unshared_gene_data_results.append(phage2_all_gene_count)
    shared_unshared_gene_data_results.append(phage2_shared_gene_count)    
    shared_unshared_gene_data_results.append(phage2_unshared_gene_count)
    shared_unshared_gene_data_results.append(round(phage2_all_genes_ave_length,4))
    shared_unshared_gene_data_results.append(round(phage2_shared_genes_ave_length,4))
    shared_unshared_gene_data_results.append(round(phage2_unshared_genes_ave_length,4))    
    shared_unshared_gene_data_results.append(phage2_all_genes_total_length)
    shared_unshared_gene_data_results.append(phage2_shared_genes_total_length)
    shared_unshared_gene_data_results.append(phage2_unshared_genes_total_length)
    shared_unshared_gene_data_results.append(round(phage2_all_gene_GC_content,4))
    shared_unshared_gene_data_results.append(round(phage2_shared_gene_GC_content,4))    
    shared_unshared_gene_data_results.append(round(phage2_unshared_gene_GC_content,4))

    output_writer.writerow(shared_unshared_gene_data_results)



    
    

    
    

    
    











output_csvfile.close()



print("Shared/unshared gene analysis completed.")



