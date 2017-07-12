#To retrieve ncbi taxonomy with a list of taxon names
#Travis Mavrich
#20160925



#Import modules
import os, sys, time, csv
from ete3 import NCBITaxa


# Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:

    input_file = sys.argv[1]
    
except:

    print("\n\n\
        This is a Python3 script to retrieve taxonomy data from the ncbi database\n\
        Execute script in any working directory\n\
        Script requires one argument(s):\n\n\
        \n\
        \n\
        First argument: file of taxon names with the following columns (csv-formatted)\n\
        1 = taxon name\
        \n")

    sys.exit(1)


#Set up the taxonomy database
print("\n\nSetting up NCBI taxonomy database")
ncbi = NCBITaxa()

#This is pre-set list of possible virus taxonomy field info that is listed as 'no rank' but which
#reflects virus DNA type
virus_type_dict = {"dsDNA viruses, no RNA stage":"dsDNA","ssDNA viruses":"ssDNA","ssRNA viruses":"ssRNA","dsRNA viruses":"dsRNA"}







#Create output files containing all phage header information and all feature information
date = time.strftime("%Y%m%d")


#Retrieve contents of the input file
input_handle = open(input_file,"r")
input_reader = csv.reader(input_handle)
taxon_set = set()
for row in input_reader:
    taxon_set.add(row[0])
input_handle.close()








#Create list of phage header names
#0 = The name was provided to retrieve data for
#1 = Full taxonomy lineage ranks, unparsed
#2 = Full taxonomy lineage, unparsed
#3 = Superkingdom
#4 = Viral type
#5 = Phylum
#6 = Class
#7 = Order
#8 = Family
#9 = Genus
#10 = Species
header_info = [\
    "Retrieval_taxon",\
    "Taxonomy_lineage ranks",\
    "Taxonomy_lineage",\
    "Superkingdom",\
    "Viral_type",\
    "Phylum",\
    "Class",\
    "Order",\
    "Family",\
    "Genus",\
    "Species"]

output_file = '%s_ncbi_taxonomy_info.csv' % date
output_handle = open(output_file,"w")
output_writer = csv.writer(output_handle)
output_writer.writerow(header_info)

error_message1 = "Unable to be retrieved"
error_message2 = "Multiple taxa with identical name"

#Convert taxon names to taxid dictionary
#Key = taxon. Value = taxid
taxid_dict = ncbi.get_name_translator(list(taxon_set))

#Iterate through the original set of taxon names. This will allow identification of taxa with no retrieved taxid
for input_taxon in taxon_set:

    print("Retrieving data for %s" % input_taxon)
    output_list = []
    taxonomy_rank_list = []
    taxonomy_list = []


    #See if the taxon was present in the ncbi database
    if input_taxon not in taxid_dict or input_taxon == "Unspecified":
        input_taxid = error_message1
        taxonomy_rank_list.append(error_message1)
        taxonomy_list.append(error_message1)
        tax_superkingdom = error_message1
        tax_viral_type = error_message1
        tax_phylum = error_message1
        tax_class = error_message1
        tax_order = error_message1
        tax_family = error_message1
        tax_genus = error_message1
        tax_species = error_message1
    

    else:
        
        input_taxid_list = taxid_dict[input_taxon]

        #If the input taxon name results in multiple taxa, then unable to parse output
        if len(input_taxid_list) > 1:

            input_taxid = error_message2
            taxonomy_rank_list.append(error_message2)
            taxonomy_list.append(error_message2)
            tax_superkingdom = error_message2
            tax_viral_type = error_message2
            tax_phylum = error_message2
            tax_class = error_message2
            tax_order = error_message2
            tax_family = error_message2
            tax_genus = error_message2
            tax_species = error_message2

        else:
            input_taxid = input_taxid_list[0]
            print(input_taxid)
            tax_superkingdom = "Unspecified"
            tax_viral_type = "Unspecified"
            tax_phylum = "Unspecified"
            tax_class = "Unspecified"
            tax_order = "Unspecified"
            tax_family = "Unspecified"
            tax_genus = "Unspecified"
            tax_species = "Unspecified"

            lineage_list = ncbi.get_lineage(input_taxid)
            lineage_rank_dict = ncbi.get_rank(lineage_list)
            taxonomy_name_dict = ncbi.get_taxid_translator(lineage_list)
        
            #Create a full record of the taxonomy rank and name data for output
            for input_taxid in lineage_list:
                taxonomy_rank_list.append(lineage_rank_dict[input_taxid])
                taxonomy_list.append(taxonomy_name_dict[input_taxid])

            #Parse specific rankings
            for output_taxid in lineage_rank_dict:
    
                #Retrieve the rank and name of the taxon
                rank = lineage_rank_dict[output_taxid]
                name_lookup = ncbi.get_taxid_translator([output_taxid])
                name = name_lookup[output_taxid]

                #If the rank matches the specified designation, save it in a variable for output
                if rank == "superkingdom":
                    tax_superkingdom = name

                elif rank == "no rank":
            
                    if name in virus_type_dict.keys():
                        tax_viral_type = virus_type_dict[name]

                elif rank == "phylum":
                    tax_phylum = name

                elif rank == "class":
                    tax_class = name

                elif rank == "order":
                    tax_order = name

                elif rank == "family":
                    tax_family = name

                elif rank == "genus":
                    tax_genus = name

                elif rank == "species":
                    tax_species = name

                else:
                    pass

    #Append output list
    output_list.append(input_taxon)
    output_list.append(str(taxonomy_rank_list))
    output_list.append(str(taxonomy_list))
    output_list.append(tax_superkingdom)
    output_list.append(tax_viral_type)
    output_list.append(tax_phylum)
    output_list.append(tax_class)
    output_list.append(tax_order)
    output_list.append(tax_family)
    output_list.append(tax_genus)
    output_list.append(tax_species)
    output_writer.writerow(output_list)


#Close script.
output_handle.close()
print("\n\n\nTaxonomy retrieval script completed.")
