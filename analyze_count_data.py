#Python3 script to analyze horizontal gene transfer data from Count software
#Travis Mavrich
#20170320



#Import built-in modules
import sys, os, csv

#Import third-party modules
from ete3 import Tree


#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    tree_file = sys.argv[1]
    hgt_file = sys.argv[2]

except:


    print("\n\n\
        This is a Python3 script to analyze horizontal gene transfer data from Count software\n\
        Execute script in any working directory\n\
        Script requires two input files:\n\n\
        First file - tree file (tab delimited):\n\
            Newick format\n\
            All internal nodes must be labeled and match nodes in count output\n\n\
        Second file - HGT data file (csv-formatted):\n\
            0 Node name in tree\n\
            1 Individual pham gains\n\
            2 Individual pham losses\n\
            \n\n\
        The following file(s) are generated:\n\
        \n\
        1. Cumulative pham gain and loss between tree leaves (csv-formatted):\n\
            0 Phage1\n\
            1 Phage2\n\
            2 Cumulative branch distance\n\
            3 Cumulative pham gains\n\
            4 Cumulative pham losses\n\
        ")
        
    sys.exit(1)
    
    
    




    
    
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






#Verify the hgt data file exists

#Expand the path if it references the home directory
if hgt_file[0] == "~":
    hgt_file = home_dir + hgt_file[1:]

if hgt_file[0] != "/":
    hgt_file = working_dir + '/' + hgt_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
hgt_file = os.path.abspath(hgt_file)


if os.path.exists(hgt_file) == False:
    print("\n\nInvalid HGT data file path.\n\n")
    sys.exit(1)

file_name = input("Count analysis output file name: ")

#Import tree data
tree = Tree(tree_file,format=3)





#Import name data into dictionary
#Key = pham
#Value = list of gains and losses
file_object = open(hgt_file,"r")
file_reader = csv.reader(file_object)
hgt_data_dict = {}
for row in file_reader:

    row[1] = int(row[1])
    row[2] = int(row[2])
    hgt_data_dict[row[0]] = row
file_object.close()


#Compute pairwise gains and losses
pairwise_comparison_data_list = []
for leaf1 in tree:

    for leaf2 in tree:
    
        #If both leaves are the same, then skip
        if leaf1.name != leaf2.name:

            print("\n\n\nCollecting data for leaves: %s, %s" %(leaf1.name,leaf2.name))
            leaf1_leaf2_common_ancestor = tree.get_common_ancestor(leaf1.name,leaf2.name)
            print("Common ancestor: " + leaf1_leaf2_common_ancestor.name)

            leaf1_leaf2_dist_sum = 0
            leaf1_leaf2_gain_sum = 0
            leaf1_leaf2_loss_sum = 0

            #For both leaves, trace back to the common ancestor and collect data
            for leaf in [leaf1,leaf2]:
                active_node = leaf
                print("\nTracing back leaf node: %s" %leaf.name)
                print("Active node: " + active_node.name)

                while active_node.name != leaf1_leaf2_common_ancestor.name:

                    active_node_hgt_data = hgt_data_dict[active_node.name]

                    print("Active node hgt gain: " + str(active_node_hgt_data[1]))
                    print("Active node hgt loss: " + str(active_node_hgt_data[2]))
                    print("Active node distance: " + str(active_node.dist))

                    leaf1_leaf2_dist_sum += active_node.dist
                    leaf1_leaf2_gain_sum += active_node_hgt_data[1]
                    leaf1_leaf2_loss_sum += active_node_hgt_data[2]



                    active_node = active_node.up
                    print("New active node: " + active_node.name)

            #Now that all nodes for both leaves have been visited and data has been collected, store it in master list
            data_list = [leaf1.name,\
                            leaf2.name,\
                            round(leaf1_leaf2_dist_sum,4),\
                            leaf1_leaf2_gain_sum,\
                            leaf1_leaf2_loss_sum]

            print(data_list)    
            pairwise_comparison_data_list.append(data_list)


#Export pairwise data
print("Exporting pairwise data...")
columns = ['phage1',\
            'phage2',\
            'distance',\
            'pham_gain_sum',\
            'pham_loss_sum']


csvfile = open(file_name,'w')
writer = csv.writer(csvfile)
writer.writerow(columns)
for line in pairwise_comparison_data_list:
    writer.writerow(line)
csvfile.close()



#End script
print("HGT analyss script completed.")



















