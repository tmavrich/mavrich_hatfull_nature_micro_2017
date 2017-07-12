#Python3 script to analyze Viral Orthologous Group (VOG) data
#Travis Mavrich



#Import built-in modules
import time, sys, os, csv, statistics







#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    accession_file = sys.argv[1]
    vog_file = sys.argv[2]

except:


    print("\n\n\
        This is a Python3 script to analyze viral orthologous genes (VOG) data \n\
        Execute script in any working directory\n\
        Script requires two input files:\n\n\
        First file: accession file (csv-formatted)\n\
            0 Phage identifier\n\
            1 Accession\n\
            2 Total number of genes\n\
            3 Total number of VOGs\n\n\
            \n\
        Second file: VOG file (tab-delimited)\n\
            0 Unparsed VOG data\n\
            \n\
            \n")




    sys.exit(1)



#Expand home and working directory
home_dir = os.path.expanduser('~')
working_dir = os.path.abspath('.')




#Verify the accession data file exists

#Expand the path if it references the home directory
if accession_file[0] == "~":
    accession_file = home_dir + accession_file[1:]

if accession_file[0] != "/":
    accession_file = working_dir + '/' + accession_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
accession_file = os.path.abspath(accession_file)


if os.path.exists(accession_file) == False:
    print("\n\nInvalid genome data file path.\n\n")
    sys.exit(1)





#Verify the vog data file exists

#Expand the path if it references the home directory
if vog_file[0] == "~":
    vog_file = home_dir + vog_file[1:]

if vog_file[0] != "/":
    vog_file = working_dir + '/' + vog_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
vog_file = os.path.abspath(vog_file)

if os.path.exists(vog_file) == False:
    print("\n\nInvalid VOG data file path.\n\n")
    sys.exit(1)




date = time.strftime("%Y%m%d")



#Import raw data
print("Importing data")



#Open accession data
accession_file_handle = open(accession_file,"r")
accession_file_reader = csv.reader(accession_file_handle)
genome_data_list = []
for row in accession_file_reader:
    row[2] = int(row[2])
    row[3] = int(row[3])
    genome_data_list.append(row)
accession_file_handle.close()



#Open VOG data
vog_file_handle = open(vog_file,"r")
vog_file_reader = csv.reader(vog_file_handle,delimiter='\t')
vog_data_list = []
for row in vog_file_reader:

    split_data = row[0].split(":")
    vog = split_data[0]
    accession = split_data[1]

    vog_data_list.append([vog,accession])
vog_file_handle.close()









#Create dictionaries, sets, and lists to perform computations


#Create genome data dictionary
#Key = accession
#Value = list of genome data
genome_data_dict = {}
for row in genome_data_list:

    if row[1] in genome_data_dict.keys():
        print("Error: %s is a duplicate accession. Unable to proceed." %row[1])
        sys.exit(1)

    genome_data_dict[row[1]] = row



#Create accession-VOG dictionary
#Key = accession
#Value = list of VOGs
accession_vog_dict = {}
for accession in genome_data_dict.keys():


    vog_list = []
    for row in vog_data_list:
        if accession == row[1]:
            vog_list.append(row[0])
    accession_vog_dict[accession] = vog_list




#Verify VOG data matches genome data
for accession in genome_data_dict.keys():


    genome_data_list = genome_data_dict[accession]
    vog_sum1 = genome_data_list[3]
    matched_vog_list = accession_vog_dict[accession]
    matched_vog_set = set(matched_vog_list)

    if vog_sum1 != len(matched_vog_set):
        print("Error for accession %s: VOG sum of %s in genome data does not match VOG sum of %s in VOG data."\
         %(accession,vog_sum1,len(matched_vog_set)))
        input("Press ENTER to continue...")
    else:
        print("VOG sums match for accession %s." %accession)











#Pairwise shared VOG proportions
print("Computing pairwise shared VOG proportions")


#This list will stored the proportions of VOGs shared between phages
proportion_data_list = []


error_list = []
for phage1 in genome_data_dict.keys():


    phage1_genome_data = genome_data_dict[phage1]
    phage1_name = phage1_genome_data[0]
    phage1_accession = phage1_genome_data[1]
    phage1_total_genes = phage1_genome_data[2]
    phage1_total_vogs = phage1_genome_data[3]
    phage1_vog_list = accession_vog_dict[phage1]
    phage1_vog_set = set(phage1_vog_list)




    #Verify the gene tallies add up
    if phage1_total_genes < len(phage1_vog_list):
        error_list.append("Error for accession %s: %s total genes %s does not match total VOG genes %s." %(phage1_accession,phage1_name,phage1_total_genes,len(phage1_vog_list)))

    for phage2 in genome_data_dict.keys():


        #If both phages to compare are the same, then skip this genome
        if phage1 == phage2:
            print("Phages %s and %s are the same." % (phage1,phage2))
            continue

        #Reset this variable for every pairwise comparison
        shared_unshared_gene_data_results = []

        phage2_genome_data = genome_data_dict[phage2]
        phage2_name = phage2_genome_data[0]
        phage2_accession = phage2_genome_data[1]
        phage2_total_genes = phage2_genome_data[2]
        phage2_total_vogs = phage2_genome_data[3]
        phage2_vog_list = accession_vog_dict[phage2]
        phage2_vog_set = set(phage2_vog_list)

        intersection_set = phage1_vog_set & phage2_vog_set
        phage1_unshared_vog_set = phage1_vog_set - intersection_set
        phage2_unshared_vog_set = phage2_vog_set - intersection_set


        #Since not all genes are given a gene VOG, it is not possible to
        #compute the exact equivalent values as with pham data.
        #As an alternative, compute the total number of genes in
        #each genome that have been grouped into any of the VOGs that are
        #shared between the two genomes.
        #Now that total gene numbers are used, and not VOG numbers,
        #gene dissimilarity is computed using the total number of genes in the genome
        
        phage1_shared_vog_gene_tally = 0
        phage1_unshared_vog_gene_tally = 0
        for vog in phage1_vog_list:
            if vog in intersection_set:
                phage1_shared_vog_gene_tally += 1
            else:
                phage1_unshared_vog_gene_tally += 1


        phage2_shared_vog_gene_tally = 0
        phage2_unshared_vog_gene_tally = 0
        for vog in phage2_vog_list:
            if vog in intersection_set:
                phage2_shared_vog_gene_tally += 1
            else:
                phage2_unshared_vog_gene_tally += 1


        phage1_unshared_other_gene_tally = phage1_total_genes - phage1_shared_vog_gene_tally - phage1_unshared_vog_gene_tally
        phage2_unshared_other_gene_tally = phage2_total_genes - phage2_shared_vog_gene_tally - phage2_unshared_vog_gene_tally

        phage1_shared_vog_gene_proportion = phage1_shared_vog_gene_tally/phage1_total_genes
        phage2_shared_vog_gene_proportion = phage2_shared_vog_gene_tally/phage2_total_genes
        ave_proportion = (phage1_shared_vog_gene_proportion + phage2_shared_vog_gene_proportion)/2
        gene_content_dissimilarity = 1 - ave_proportion


        proportion_data_list.append([\
                                    phage1_name,\
                                    phage1_accession,\
                                    phage1_total_genes,\
                                    phage1_total_vogs,\
                                    len(phage1_unshared_vog_set),\
                                    phage1_shared_vog_gene_tally,\
                                    phage1_unshared_vog_gene_tally,\
                                    phage1_unshared_other_gene_tally,\
                                    round(phage1_shared_vog_gene_proportion,2),\
                                    \
                                    phage2_name,\
                                    phage2_accession,\
                                    phage2_total_genes,\
                                    phage2_total_vogs,\
                                    len(phage2_unshared_vog_set),\
                                    phage2_shared_vog_gene_tally,\
                                    phage2_unshared_vog_gene_tally,\
                                    phage2_unshared_other_gene_tally,\
                                    round(phage2_shared_vog_gene_proportion,2),\
                                    \
                                    len(intersection_set),\
                                    round(ave_proportion,2),\
                                    round(gene_content_dissimilarity,2)])




#Output shared VOG proportion data
print("Exporting proportion data to csv file")
columns = [\
            'phage1_name',\
            'phage1_accession',\
            'phage1_number_of_genes',\
            'phage1_number_of_vogs',\
            'phage1_number_of_unshared_vogs',\
            'phage1_number_of_shared_vog_genes',\
            'phage1_number_of_unshared_vog_genes',\
            'phage1_number_of_unshared_other_genes',\
            'phage1_shared_vog_gene_proportion',\
            \
            'phage2_name',\
            'phage2_accession',\
            'phage2_number_of_genes',\
            'phage2_number_of_vogs',\
            'phage2_number_of_unshared_vogs',\
            'phage2_number_of_shared_vog_genes',\
            'phage2_number_of_unshared_vog_genes',\
            'phage2_number_of_unshared_other_genes',\
            'phage2_shared_vog_gene_proportion',\
            \
            'number of shared vogs',\
            'average_shared_vog_gene_proportion',\
            'gene_content_dissimilarity']


output_file_name =  '%s_shared_vog_proportion_data.csv' %date
output_file_handle = open(output_file_name,'w')
output_file_writer = csv.writer(output_file_handle)
output_file_writer.writerow(columns)

for line in proportion_data_list:
    output_file_writer.writerow(line)
output_file_handle.close()




#Print out list of phages with errors
if len(error_list) > 0:
    print("Error: some phages have inconsistent gene tallies:")
    for element in error_list:
        print(element)

print("\n\n\nVOG analysis completed.")
