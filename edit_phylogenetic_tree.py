#Python3 script to manipulate phylogenetic tree data
#Travis Mavrich



#Import built-in modules
import time, sys, os, csv

#Import third-party modules
from ete3 import Tree


#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    tree_file = sys.argv[1]
    name_file = sys.argv[2]

except:


    print("\n\n\
        This is a Python3 script to manipulate phylogenetic tree data\n\
        Execute script in any working directory\n\
        Script requires two input files:\n\n\
        First file: tree file (tab delimited)\n\
            \n\
        Second file: name file (csv-formatted):\n\
            0 Taxon name in tree\n\
            1 New taxon name\n\
            (Note: if you don't want to convert taxa names, enter 'none'.\
            \n\
            \n")



    sys.exit(1)



#Determine what the user wants to do
print("\n\nChoose one option from below:")
print("1 = Do you want to rename tree leaf names?")
print("2 = Do you want to label internal nodes?")
print("3 = Do you want to parse elements of the tree?")

response_valid = False
response = input()

while response_valid == False:
    if response == "1" or \
        response == "2" or \
        response == "3":

        response = int(response)
        response_valid = True

    else:
        print("Invalid response.")




#Determine what the format of the input tree is in
#Format 0 = flexible, allows missing information
#Format 2 = all branch names, all leaf names, internal branch support
#Format 3 = all branch names, all leaf names
#If you import a tree that contains node names in the place of supports, and don't specify the format, it throws an error.
#Default support value = 1.


print("\n\nChoose the input tree format:")
print("0 = It looks for support values, but doesn't need them. (Choose this if unsure)")
print("2 = branch lengths, leaf node names, and branch supports")
print("3 = branch lengths, leaf node names, internal node names, no branch supports")
print("5 = branch lengths, leaf node names")
print("8 = leaf node names, internal node names")
print("9 = leaf node names")

input_format_valid = False

while input_format_valid == False:

    input_format = input()


    if input_format == "0" or \
        input_format == "2" or \
        input_format == "3" or \
        input_format == "5" or \
        input_format == "8" or \
        input_format == "9":

        input_format = int(input_format)
        input_format_valid = True

    else:
        print("Invalid response.")




#Expand home and working directory
home_dir = os.path.expanduser('~')
working_dir = os.path.abspath('.')


#Verify the tree data file exists

#Expand the path if it references the home directory
tree_basename = tree_file.split('.')[0]


if tree_file[0] == "~":
    tree_file = home_dir + tree_file[1:]

if tree_file[0] != "/":
    tree_file = working_dir + '/' + tree_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
tree_file = os.path.abspath(tree_file)


if os.path.exists(tree_file) == False:
    print("\n\nInvalid tree data file path.\n\n")
    sys.exit(1)











#Import tree data
try:
    tree = Tree(tree_file,format=input_format)
except:
    print("Error: unable to import tree with the specified format.")
    print("Exiting script.")
    sys.exit(1)


#Rename tree leaf names if selected by user
if response == 1:

    if name_file.lower() != "none":
        #Expand the path if it references the home directory
        if name_file[0] == "~":
            name_file = home_dir + name_file[1:]

        if name_file[0] != "/":
            name_file = working_dir + '/' + name_file

        #Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
        name_file = os.path.abspath(name_file)

        if os.path.exists(name_file) == False:
            print("\n\nInvalid name data file path.\n\n")
            sys.exit(1)
    else:
        print("Error: no name file was provided.")
        print("Exiting script.")
        sys.exit(1)


    #Import name data
    file_object = open(name_file,"r")
    file_reader = csv.reader(file_object)
    name_conversion_dict = {}
    for row in file_reader:
        name_conversion_dict[row[0]] = row[1]
    file_object.close()


    #Replace all leaf names
    for leaf in tree:

        try:
            leaf.name = name_conversion_dict[leaf.name]
        except:
            input("Error: unable to replace %s." %leaf.name)




#Label internal nodes if selected by user
if response == 2:

    node_label_num = 1
    for node in tree.traverse("postorder"):

        if node.is_leaf() == False and node.name == "":
            node.name = "node%s" %node_label_num
            node_label_num += 1
            print("Internal node changed to %s" %node.name)
        else:
            print("nothing changed")
            print(node.name)

        if node.is_root() == True:
            print("Root is: %s" %node.name)


#Parse tree elements if selected by user
if response == 3:

    node_data_list = []
    for node in tree.traverse("postorder"):
        node_data_list.append([node.name,node.dist,node.support])


    #Export node data
    print("Exporting parsed tree data...")
    columns = ['node_name',\
                'branch_length',\
                'branch_support']


    csvfile = open('parsed_tree_data.csv','w')
    writer = csv.writer(csvfile)
    writer.writerow(columns)
    for line in node_data_list:
        writer.writerow(line)
    csvfile.close()







#Export the renamed tree based on selections
#If you specify format=3, then support values are lost, and node names are
#outputted in the support position, even if there are no node names (default = "NoName")
#If you don't specify format, it looks to output anything available.
#Not sure if it gives preference to support or node names


if response == 1 or response == 2:


    print("\n\nChoose the output tree format:")
    print("0 = It looks for support values, but doesn't need them. (Choose this if unsure)")
    print("2 = branch lengths, leaf node names, and branch supports")
    print("3 = branch lengths, leaf node names, internal node names, no branch supports")
    print("5 = branch lengths, leaf node names")
    print("8 = leaf node names, internal node names")
    print("9 = leaf node names")

    output_format_valid = False

    while output_format_valid == False:

        output_format = input()

        if output_format == "0" or \
            output_format == "2" or \
            output_format == "3" or \
            output_format == "5" or \
            output_format == "8" or \
            output_format == "9":

            output_format = int(output_format)
            output_format_valid = True

        else:
            print("Invalid response.")


    #Format_root_node will output the root node, if it is present.
    #Otherwise, it doesn't seem to care if the root is absent.
    tree.write(outfile=tree_basename + "_modified.txt",format=output_format,format_root_node=True)





#End script
print("\n\nTree manipulation script completed.")
