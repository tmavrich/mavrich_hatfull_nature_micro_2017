#Python3 script to process Mash output
#It removes all duplicate reactions and trims the names





#Import modules
import csv, os, sys


def trim(phage_name):
    phage_name = phage_name.split("/")[-1]
    phage_name = phage_name[:-6]
    return(phage_name)



#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    input_file = sys.argv[1]
    output_file = sys.argv[2]

except:


    print("\n\n\
        This is a Python3 script to process Mash output files, in preparation for downstream analyses\n\
        Script requires two arguments:\n\n\
        Execute script in any working directory\n\
        First argument: input file, which is the Mash analysis output (tab-delimited)\n\
            0 = reference_path\n\
            1 = query_path\n\
            2 = mash distance\n\
            3 = mash p-value\n\
            4 = kmer count\n\
            \n\
        Second argument: filename to store the output data (csv-formatted)\n\
            0 = reference\n\
            1 = query\n\
            2 = mash distance\n\
            3 = mash p-value\n\
            4 = kmer count\n\
            5 = ref_query\n\
            \n")
    sys.exit(1)


#Open the mash output file
#0 = reference_path
#1 = query_path
#2 = mash distance
#3 = mash p-value
#4 = kmer count
mash_data_list = []
mash_data_handle = open(input_file, 'r')
mash_reader = csv.reader(mash_data_handle,delimiter='\t')
for row in mash_reader:
    mash_data_list.append(row)

mash_data_handle.close



#Parse every line and remove if it contains duplicate data
mash_dict = {}
total = 0
added = 0
not_added = 0
same_phage = 0
print("Comparing lines of data...")
for line in mash_data_list:
    total += 1
    
    #Trim off the preceding path text and ".fasta" suffixes.
    reference = trim(line[0])
    query = trim(line[1])
    
    
    #Remove if reference and query are the same.
    if reference == query:
        print("Reference and query are the same: %s, %s" % (reference,query))
        same_phage += 1
        continue
    
    #Create a reference_query name combination
    ref_query_list = [reference,query]
    ref_query_list.sort()
    ref_query = ref_query_list[0] + "_" + ref_query_list[1]
    print(ref_query)
    
    #If the reference_query name combination is not already in the dictionary, add it.
    #Otherwise, verify that the the matching reference_query name combination contains matching data.
    #(which is should). If it doesn't match, throw an error.
    if ref_query not in mash_dict:
        line[0] = reference
        line[1] = query
        line.extend([ref_query])
        mash_dict[ref_query] = line
        added += 1
        print(line)
    else:
        not_added += 1
        compare_line = mash_dict[ref_query]
        if (compare_line[2] != line[2] or compare_line[3] != line[3] or compare_line[4] != line[4]):
            print("Error: non-matching data rows")
            print(line)
            print(compare_line)
        else:
            pass



print("Total input rows: %s" % total)
print("Added input rows: %s" % added)
print("Not added: %s" % not_added)
print("Same phage: %s" % same_phage)

#Output into a new csv file
print("Outputting results to file...")
output_file_handle = open(output_file,"w")
writer = csv.writer(output_file_handle)
column_headers = ["reference","query","distance","p-value","kmer_count","ref_query"]
writer.writerow(column_headers)
for key in mash_dict:
    writer.writerow(mash_dict[key])
output_file_handle.close()    
    
    
    
    
    
    
    

