#To extract all features and/or header information in a genbank file
#Travis Mavrich


#Import modules
import os, sys, time, csv
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from ete3 import NCBITaxa



# Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:

    genbank_files_dir = sys.argv[1]

except:

    print("\n\n\
        This is a Python3 script to extract feature information from Genbank-formatted files\n\
        Execute script in any working directory\n\
        Script requires one argument(s):\n\n\
        \n\
        \n\
        First argument: directory of genomes (Genbank-formatted)\n\
        \n\n\n\n\
        The script outputs three files:\n\
        1. Summary of genome information (csv-formatted)\n\
        2. Summary of CDS and tRNA feature information (csv-formmated)\n\
        3. Concatenated list of gene translations per genome (fasta-formatted)\n\
         \n")

    sys.exit(1)


#Set up the taxonomy database
print("\n\nSetting up NCBI taxonomy database")
ncbi = NCBITaxa()

#This is pre-set list of possible virus taxonomy field info that is listed as 'no rank' but which
#reflects virus DNA type
virus_type_dict = {"dsDNA viruses, no RNA stage":"dsDNA","ssDNA viruses":"ssDNA","ssRNA viruses":"ssRNA","dsRNA viruses":"dsRNA"}


#Set up protein alphabet to verify sequence integrity
protein_alphabet_set = set(IUPAC.ExtendedIUPACProtein.letters)


#Create list of phage header names
#0 = Filename
#1 = Phage_Name = phageName
#2 = Record_Name
#3 = Accession_Number = accessionNum
#4 = Record_ID
#5 = Record_Definition
#6 = Record_Source
#7 = Record_Organism
#8 = Source_Feature_Organism
#9 = Source_Feature_Host
#10 = Source_Feature_Lab_Host
#11 = RecordDef_Genus
#12 = RecordSource_Genus
#13 = RecordOrganism_Genus
#14 = SourceFeatureOrganism_Genus
#15 = SourceFeatureHost_Genus
#16 = SourceFeatureLabHost_Genus
#17 = GenomeSize = seqLength
#18 = Genus_comparison
#19 = Full taxonomy lineage
#20 = Superkingdom
#21 = Viral type
#22 = Phylum
#23 = Class
#24 = Order
#25 = Family
#26 = Genus
#27 = Species
#28 = Species verification
header_info = [\
    "Filename",\
    "Phage_Name",\
    "Record_Name",\
    "AccessionNumber",\
    "Record_ID",\
    "Record_Defintion",\
    "Record_Source",\
    "Record_Organism",\
    "Source_Feature_Organism",\
    "Source_Feature_Host",\
    "Source_Feature_Lab_Host",\
    "RecordDef_Genus",\
    "RecordSource_Genus",\
    "RecordOrganism_Genus",\
    "SourceFeatureOrganism_Genus",\
    "SourceFeatureHost_Genus",\
    "SourceFeatureLabHost_Genus",\
    "GenomeSize",\
    "GenusComparison",\
    "Taxonomy_lineage",\
    "Superkingdom",\
    "Viral_type",\
    "Phylum",\
    "Class",\
    "Order",\
    "Family",\
    "Genus",\
    "Species",\
    "Species_verification"]


#Create list of feature header names
#0 = basename
#1 = phageName
#2 = locus_tag
#3 = gene number
#4 = startCoord (uncorrected)
#5 = stopCoord
#6 = typeID
#7 = translation table
#8 = translation
#9 = orientation
#10 = feature_product
#11 = feature_function
#12 = feature_note
feature_info = [\
"Filename",\
"Phage_Name",\
"Locus_Tag",\
"Gene_Number",\
"StartCoord (uncorrected)",\
"StopCoord ",\
"TypeID",\
"Translation_table",\
"Translation",\
"Orientation",\
"Product",\
"Function",\
"Note"]




#Create output files containing all phage header information and all feature information
date = time.strftime("%Y%m%d")


phage_header_file = '%s_phage_header_info.csv' % date
phage_header_handle = open(phage_header_file,"w")
phage_header_writer = csv.writer(phage_header_handle)
phage_header_writer.writerow(header_info)


phage_feature_file = '%s_phage_feature_info.csv' % date
phage_feature_handle = open(phage_feature_file,"w")
phage_feature_writer = csv.writer(phage_feature_handle)
phage_feature_writer.writerow(feature_info)

concatenated_translations_file = '%s_concatenated_translations.csv' % date
concatenated_translations_handle = open(concatenated_translations_file,"w")
concatenated_translations_writer = csv.writer(concatenated_translations_handle)




#Iterate through each Genbank file in the directory
for filename in os.listdir(genbank_files_dir):

    print(filename)
    seq_record = SeqIO.read(os.path.join(genbank_files_dir,filename), 'genbank')


    #Initialize variables.
    #This ensures no data is carried over between files for some variables that may not be present in all files
    basename = filename.split('.')[0]
    feature_source_organism = ""
    feature_source_host = ""
    feature_source_lab_host = ""
    phage_data_list = []
    all_features_data_list = []
    translation_list = []
    feature_source_list = []
    total_protein_sequence_length = 0
    taxid = ""
    tax_superkingdom = "Unable to be retrieved"
    tax_viral_type = "Unable to be retrieved"
    tax_phylum = "Unable to be retrieved"
    tax_class = "Unable to be retrieved"
    tax_order = "Unable to be retrieved"
    tax_family = "Unable to be retrieved"
    tax_genus = "Unable to be retrieved"
    tax_species = "Unable to be retrieved"
    species_check = "Unable to verify"




    #Retrieve the header information


    #Record organism, Phage Name
    try:
        record_organism = seq_record.annotations["organism"]
        phageName = record_organism.split(' ')[-1]
    except:
        record_organism = ""
        phageName = ""


    record_organism_genus = record_organism.split(' ')[0]


    #Sequence, length
    try:
        seqLength = len(seq_record.seq)

    except:
        seqLength = ""

    #File header fields are retrieved to be able to check phageName and HostStrain typos
    #The Accession field, with the appended version number, is stored as the record.id
    #The Locus name at the top of the file is stored as the record.name
    #The base accession number, without the version, is stored in the 'accession' annotation list

    try:
        record_name = str(seq_record.name)
    except:
        record_name = ""

    try:
        record_id = str(seq_record.id)
    except:
        record_id = ""

    try:
        record_def = str(seq_record.description)
    except:
        record_def = ""

    record_def_genus = record_def.split(' ')[0]


    try:
        record_source = str(seq_record.annotations["source"])
    except:
        record_source = ""

    record_source_genus = record_source.split(' ')[0]

    try:
        accessionNum = seq_record.annotations["accessions"][0]
    except:
        accessionNum = ""

    try:
        taxonomy_list = seq_record.annotations["taxonomy"]
    except:
        taxonomy_list = [""]



    #Parse info from the organism's features
    for feature in seq_record.features:

        #Some genbank files have multiple source features. This complicates taxonomy parsing.
        #Collect all source features and process them afterwards.
        if feature.type == "source":
            feature_source_list.append(feature)

        #Only parse if it is a CDS or tRNA feature
        elif feature.type == "CDS" or feature.type == "tRNA":

            typeID = feature.type

            #This will store all data for this feature that will be imported
            feature_data_list = []


            #Feature_locus_tag
            try:
                feature_locus_tag = feature.qualifiers["locus_tag"][0]
            except:
                feature_locus_tag = ""

            #Gene number
            try:
                feature_gene_num = feature.qualifiers['gene'][0]
            except:
                feature_gene_num = ""

            #Orientation
            if feature.strand == 1:
                orientation = "Forward"
            elif feature.strand == -1:
                orientation = "Reverse"
            #ssRNA phages
            elif feature.strand is None:
                orientation = "None"
            else:
                orientation = ""

            #Gene boundary coordinates
            #Compound features are tricky to parse.
            if str(feature.location)[:4] == "join":

                #Skip this compound feature if it is comprised of more than two features (too tricky to parse).
                if len(feature.location.parts) > 2:

                    strStart = ""
                    strStop = ""

                else:

                    #Retrieve compound feature positions based on strand
                    if feature.strand == 1:

                        strStart = str(feature.location.parts[0].start)
                        strStop = str(feature.location.parts[1].end)

                    elif feature.strand == -1:

                        strStart = str(feature.location.parts[1].start)
                        strStop = str(feature.location.parts[0].end)

                    #If strand is None...
                    else:
                        strStart = ""
                        strStop = ""

            else:
                strStart = str(feature.location.start)
                strStop = str(feature.location.end)



            #Translation
            try:
                translation = feature.qualifiers["translation"][0].upper()

            except:
                translation = ""

            total_protein_sequence_length += len(translation)

            #Compile list of translations for downstream analysis if there are no apparent sequence errors
            amino_acid_set = set(translation)
            amino_acid_error_set = amino_acid_set - protein_alphabet_set
            if len(amino_acid_error_set) == 0:
                translation_list.append(translation)


            #Translation table used
            try:
                feature_transl_table = feature.qualifiers["transl_table"][0]
            except:
                feature_transl_table = ""


            #Gene descriptions
            try:
                feature_product = feature.qualifiers["product"][0].strip()
            except:
                feature_product = ""

            try:
                feature_function = feature.qualifiers["function"][0].strip()
            except:
                feature_function = ""

            try:
                feature_note = feature.qualifiers["note"][0].strip()
            except:
                feature_note = ""



            #Now that it has acquired all gene feature info,
            #create list of gene data and append to list of all gene feature data
            #0 = basename
            #1 = phageName
            #2 = locus_tag
            #3 = gene number
            #4 = startCoord (uncorrected)
            #5 = stopCoord
            #6 = typeID
            #7 = translation table
            #8 = translation
            #9 = orientation
            #10 = feature_product
            #11 = feature_function
            #12 = feature_note
            feature_data_list.append(basename)
            feature_data_list.append(phageName)
            feature_data_list.append(feature_locus_tag)
            feature_data_list.append(feature_gene_num)
            feature_data_list.append(strStart)
            feature_data_list.append(strStop)
            feature_data_list.append(typeID)
            feature_data_list.append(feature_transl_table)
            feature_data_list.append(translation)
            feature_data_list.append(orientation)
            feature_data_list.append(feature_product)
            feature_data_list.append(feature_function)
            feature_data_list.append(feature_note)
            all_features_data_list.append(feature_data_list)

        #For all other features, just skip
        else:
            continue



    #Try to parse source feature(s)
    for feature in feature_source_list:

        try:
            feature_source_organism = str(feature.qualifiers["organism"][0])
        except:
            feature_source_organism = ""

        #Only try to parse information about host if the source feature matches the record organism
        if feature_source_organism == record_organism:

            try:
                #The db_xref field can have multiple references. Make sure it is for the taxon id number
                db_xref_list = feature.qualifiers["db_xref"]
                for db_xref in db_xref_list:
                    db_xref_split = db_xref.split(':')
                    if db_xref_split[0] == "taxon":
                        taxid = db_xref_split[1]
            except:
                taxid = ""

            try:
                feature_source_host = str(feature.qualifiers["host"][0])
            except:
                feature_source_host = ""

            try:
                feature_source_lab_host = str(feature.qualifiers["lab_host"][0])
            except:
                feature_source_lab_host = ""




    #Now that source feature info has been retrieved, taxonomic data and genus info can be parsed
    feature_source_organism_genus = feature_source_organism.split(' ')[0]
    feature_source_host_genus = feature_source_host.split(' ')[0]
    feature_source_lab_host_genus = feature_source_lab_host.split(' ')[0]


    genus_set = set()

    if len(record_def_genus) > 0:
        genus_set.add(record_def_genus)
    if len(record_source_genus) > 0:
        genus_set.add(record_source_genus)
    if len(record_organism_genus) > 0:
        genus_set.add(record_organism_genus)
    if len(feature_source_organism_genus) > 0:
        genus_set.add(feature_source_organism_genus)
    if len(feature_source_host_genus) > 0:
        genus_set.add(feature_source_host_genus)
    if len(feature_source_lab_host_genus) > 0:
        genus_set.add(feature_source_lab_host_genus)

    #See if all genus data are the same
    if len(genus_set) == 0:
        genus_comparison = "error"
    elif len(genus_set) == 1:
        genus_comparison = record_def_genus
    else:
        genus_comparison = "multiple"


    #Retrieve taxonomic data if the file has a taxid, otherwise skip
    if taxid != "":

        #Since each organism may not have all rankings present, first set all rankings
        #to 'unspecified' as a positive indication that the taxonomy was retrieved.
        tax_superkingdom = "Unspecified"
        tax_viral_type = "Unspecified"
        tax_phylum = "Unspecified"
        tax_class = "Unspecified"
        tax_order = "Unspecified"
        tax_family = "Unspecified"
        tax_genus = "Unspecified"
        tax_species = "Unspecified"
        species_check = "Unspecified"

        lineage_list = ncbi.get_lineage(taxid)
        lineage_rank_dict = ncbi.get_rank(lineage_list)
        for taxon in lineage_rank_dict:

            #Retrieve the rank and name of the taxon
            rank = lineage_rank_dict[taxon]
            name_lookup = ncbi.get_taxid_translator([taxon])
            name = name_lookup[taxon]

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

    #As a way to double-check taxonomy, make sure species rank is the same as the organism name
    if record_organism == tax_species:
        species_check = "okay"
    else:
        print("Taxonomic species different than organism name")
        species_check = "different"


    #Append output list
    phage_data_list.append(basename)
    phage_data_list.append(phageName)
    phage_data_list.append(record_name)
    phage_data_list.append(accessionNum)
    phage_data_list.append(record_id)
    phage_data_list.append(record_def)
    phage_data_list.append(record_source)
    phage_data_list.append(record_organism)
    phage_data_list.append(feature_source_organism)
    phage_data_list.append(feature_source_host)
    phage_data_list.append(feature_source_lab_host)
    phage_data_list.append(record_def_genus)
    phage_data_list.append(record_source_genus)
    phage_data_list.append(record_organism_genus)
    phage_data_list.append(feature_source_organism_genus)
    phage_data_list.append(feature_source_host_genus)
    phage_data_list.append(feature_source_lab_host_genus)
    phage_data_list.append(seqLength)
    phage_data_list.append(genus_comparison)
    phage_data_list.append(str(taxonomy_list))
    phage_data_list.append(tax_superkingdom)
    phage_data_list.append(tax_viral_type)
    phage_data_list.append(tax_phylum)
    phage_data_list.append(tax_class)
    phage_data_list.append(tax_order)
    phage_data_list.append(tax_family)
    phage_data_list.append(tax_genus)
    phage_data_list.append(tax_species)
    phage_data_list.append(species_check)









    #Output the data to separate files...

    #Output phage header info
    phage_header_writer.writerow(phage_data_list)

    #Output feature info
    for feature_data in all_features_data_list:
        phage_feature_writer.writerow(feature_data)

    #Output translation info
    if len(translation_list) > 0:
        concatenated_translations = "".join(translation_list)
        if total_protein_sequence_length == len(concatenated_translations):
            print('.')
            concatenated_translations_writer.writerow([">" + basename])
            concatenated_translations_writer.writerow([concatenated_translations])
        else:
            print('Translation length error')




#Close script.
print("\n\n\n\nParsing script completed.")
phage_header_handle.close()
phage_feature_handle.close()
concatenated_translations_handle.close()
