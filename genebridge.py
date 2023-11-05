#!/usr/bin/env python3
#
#

###Mapping and Master list Program###
#Import librarys
import subprocess
import sys
import pkg_resources
def install(package): #Installing process for dependencies
    try:
        # Check if the package is already installed
        pkg_resources.get_distribution(package)
        print(f"{package} is already installed.")
    except pkg_resources.DistributionNotFound:
        # If the package is not installed, install it
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"{package} has been installed.")

#Import NetworkX the graphing library
install("networkx")
install("matplotlib")
install("numpy")
install("pandas")
install("scipy")
install("pyvis")
install("str2bool")
print('\n')
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from pyvis.network import Network
from str2bool import str2bool
import argparse, os, re

###USER DEFINED VARIABLES###
##################################
parser = argparse.ArgumentParser()

#fileName, edge2comm Table Input
parser.add_argument("-a", "--e2c", required=False,
                    nargs='?', default="edge2comm.txt", const="edge2comm.txt",
                    help="Edge to Community file suffix\ndefault='edge2comm.txt'")

#fileName2, comm2nodes Table Input
parser.add_argument("-b", "--c2n", required=False,
                    nargs='?', default="comm2nodes.txt", const="comm2nodes.txt",
                    help="Community to Nodes file suffix\ndefault='comm2nodes.txt'")

#outfile1, CTD Chemical-Gene Interaction Table Name
parser.add_argument("-c", "--ctd", required=False,
                    nargs='?', default="interactionsCTD", const="interactionsCTD",
                    help="CTD chemical-gene interaction file name\ndefault='interactionsCTD'")

#organismID, NCBI Taxonomy Number
parser.add_argument("-g", "--organism", required=False,
                    nargs='?', default=9606, const=9606, 
                    type=int,
                    help="organism NCBI Taxonomy ID number\ndefault=9606")

#outputHeader, Edge List (outfile3) Header Enable/Disable
parser.add_argument("-e", "--header", required=False,
                    nargs='?', default=True, const=True,
                    help="header option for the final edge list\ncalled '-e' because '-h' is help\ndeafult=True")

#chugReduction, Massive Timer Saver Enable/Disable
parser.add_argument("-f", "--nochug", required=False,
                    nargs='?', default=True, const=True,
                    help="Removes group_betweenness_centrality community measurement, stopping massive chugging\ncalled '-f' because not many letters left\ndeafult=True")

args = parser.parse_args()
##################################

#Define Input and Output Files
fileName = args.e2c
fileName2 = args.c2n
outfile1 = args.ctd #ChemicalGene Interaction Table, but more importantly, list of input chemicals
organism= args.organism #Define Taxonomy ID
outputHeader = str2bool(str(args.header)) #Toggle for having headers in the final node library
chugReduction = str2bool(str(args.nochug))
###USER DEFINED FUNCTIONS###
def findFile(x):
    directory_path = "."
    suffix_pattern = r'.*' + re.escape(x) + r'$'
    # List all files in the directory
    files_in_directory = os.listdir(directory_path)
    # Use regular expressions to filter files with the desired suffix
    matching_files = [file for file in files_in_directory if re.match(suffix_pattern, file)]
    # Print the matching file names
    for file_name in matching_files:
        print("Found matching file:", file_name)
    if len(matching_files) == 1:
        return(matching_files[0])
    else:
        print("Too many or no files were found with that suffix :,(")

fileName = findFile(fileName)
fileName2= findFile(fileName2)

##Parse data from a edgelist into a graph, including community
G = nx.read_edgelist(path=fileName, delimiter='\t', data=(('group', int),), create_using=nx.Graph)

#Dictionary for comm2node
data_dict = {}

# Open the TSV file and read it line by line
with open(fileName2, 'r') as tsvfile:
    # Iterate through each line in the TSV file
    for line in tsvfile:
        # Split the line into columns using tab as the delimiter
        columns = line.strip().split('\t')
        
        # The first column is used as the value
        value = columns[0]
        # The rest of the columns are used as keys
        keys = columns[1:]

        #skipping adding communities with less than 3 nodes.
        if len(keys) < 3:
            print('Community '+ value + ' was fewer than 3 nodes! Not added to graph nodes.')
            continue
        # Create entries in the dictionary
        for key in keys:
            #checking if the key already exist
            #if true, then keeps existing value and adds new value too
            if key in data_dict:
                data_dict[key] += ',' + value
            #else, just add that key and value
            else:
                data_dict[key] = value

nx.set_node_attributes(G,data_dict,'group')




###Calculations for master list###
#Calculate number of nodes
print("Nodes:",len(G.nodes))

#Calculate number of edges
print("Edges:",len(G.edges))

#NODE METRICS#
nM = pd.DataFrame()
#really REALLY important addition of personality to the program :)
print('\nNode Data chugging time.')

#List and dictionaries for the Statistics
N_L = list(G.nodes)
N_D = dict(G.degree)
B_C = nx.betweenness_centrality(G) #dictionary
C_C = nx.closeness_centrality(G) #dictionary
E_C = nx.eigenvector_centrality_numpy(G) #dictionary
C_Co = nx.clustering(G) #dictionary
P_R = nx.pagerank(G) #dictionary
group = nx.get_node_attributes(G,'group') #dictionary

#Adding Data into nM
nM['Node_Label'] = N_L                                      #(Example: APP) 
nM['Node_Degree'] = nM['Node_Label'].map(N_D)               #The number of interactions (edges) a gene has with other genes. It provides information about a gene's importance within the network.
nM['Betweenness_Centrality'] = nM['Node_Label'].map(B_C)    #Identifies genes that act as bridges or intermediaries between different parts of the network.
nM['Closeness_Centrality'] = nM['Node_Label'].map(C_C)      #Measures how close a gene is to all other genes in the network, indicating its potential influence.
nM['Eigenvector_Centrality'] = nM['Node_Label'].map(E_C)    #Identifies genes that are connected to other highly connected genes, indicating their importance in the network.
nM['Clustering_Coefficient'] = nM['Node_Label'].map(C_Co)   #Measures the extent to which a gene's neighbors are also connected to each other, indicating potential involvement in functional modules.
nM['PageRank'] = nM['Node_Label'].map(P_R)                  #Similar to eigenvector centrality, it quantifies the importance of a gene based on its connections.
nM['Community'] = nM['Node_Label'].map(group)               #Community(ies) a node is in. Nodes with multiple communities will be have them in chronological order seperated by a comma. (example: 2,4,27) Communities smaller than 3 will not be added!
#Outputs the file for Node Metrics
nM.to_csv('nodeMetrics.tsv', index=False, sep='\t', header= outputHeader)
#super important personality continuity
print("\nnodeMetrics.tsv saved. Begining community metrics.")


##COMMUNITY METRICS
#List and dictionaries for the Statistics FUNCTION (c2n standing for comm2node)
def cMetrics(c2n,chugReduction = True):

    cM = pd.DataFrame()
    edge_dict = nx.get_edge_attributes(G, 'group')
    z = 1 #for counting communities finished to you on chugReduction = False
    with open(c2n, 'r') as file:
        for line in file:
            # Split each line by tabs and remove leading/trailing whitespace
            values = line.strip().split('\t')
            # Extract the first value (group) and the remaining values into a list
            group = values[0]
            nodes =[value for value in values[1:]]

            # Initialize an empty list to store edge
            group_edges = []
            # Iterate through the dictionary and check for matching values
            for edge, value in edge_dict.items():
                if int(value) == int(group):
                    group_edges.append(edge)
            
            #Total Node Degree Count
            tndc_dict = dict(nx.degree(G,nodes))
            count = 0
            for degree in tndc_dict.values():
                count += int(degree)
            #Group Degree Centrality
            # d_c = nx.group_degree_centrality(G,nodes) #Doesn't get recognized occasionally, so running the function myself.
            d_c = len(set().union(*[set(G.neighbors(i)) for i in nodes]) - set(nodes))
            d_c /= len(G.nodes()) - len(nodes)

            #Group Closeness Centrality
            c_c = nx.group_closeness_centrality(G,nodes)
            if chugReduction == False:
                #Group Betweeness Centrality #VERY COMPUTER HEAVY, AS IN USE A SUPER COMPUTER#
                b_c = nx.group_betweenness_centrality(G,nodes)
                # Make a pandas dataframe for one group/community, then append that to the cM dataframe
                temp = pd.DataFrame({
                                    'Community': [group],              #(example: 1)
                                    'Node_Count': [len(nodes)],        #Community Size
                                    'Edge_Count': [len(group_edges)],  #Community Interactions
                                    'Total_Node_Degree_Count': [count],#Count of degrees in nodes
                                    'Degree_Centrality': [d_c],        #Group degree centrality of a group of nodes is the fraction of non-group members connected to group members.
                                    'Closeness_Centrality': [c_c],     #Group closeness centrality of a group of nodes is a measure of how close the group is to the other nodes in the graph.
                                    'Betweeness_Centrality': [b_c]     #Group betweenness centrality of a group of nodes is the sum of the fraction of all-pairs shortest paths that pass through any vertex in the group nodes.
                                    })
                print("That's community ",z," done.")
                z+=1
            else:
                # Make a pandas dataframe for one group/community, then append that to the cM dataframe
                temp = pd.DataFrame({
                                    'Community': [group],              #(example: 1)
                                    'Node_Count': [len(nodes)],        #Community Size
                                    'Edge_Count': [len(group_edges)],  #Community Interactions
                                    'Total_Node_Degree_Count': [count],#Count of degrees in nodes
                                    'Degree_Centrality': [d_c],        #Group degree centrality of a group of nodes is the fraction of non-group members connected to group members.
                                    'Closeness_Centrality': [c_c],     #Group closeness centrality of a group of nodes is a measure of how close the group is to the other nodes in the graph.
                                    })
            # Append the new row to the DataFrame
            cM = pd.concat([cM, temp], ignore_index=True)

    #Outputs the file for Community Metrics
    cM.to_csv('commMetrics.tsv', index=False, sep='\t', header= outputHeader)

cMetrics(fileName2,chugReduction)
print('commMetrics.tsv saved. Now onto finishing the graph visualization.')

##Ideas for adding later... shhhh don't tell anyone if you see this and it has been a few months since the last update... except me of course haha, I'll try getting back to this.##
#Take the edge2comm.txt file, copy it to main folder and rename it 'Community Edge List'
#At bottom or top, add MASTER community which accounts for all
#Add second list accounting for Shortest paths between all input chemicals, transtivity, & other network wide metrics 
##



##Dictionary for chemical labels
#chemList for input chemicals
chemList = []
of1df = pd.read_table(outfile1+'_chemical-protein.tsv') #outfile1 dataframe code
#Select for only human data(assuming human); haa stands for "I'm only Human, After All" (its a meme)
haa = of1df[of1df["OrganismID"] == organism]
##select for the column value, chemicalName, and save as a list & string
chemicalName = haa[["ChemicalName"]].drop_duplicates(keep='first')
for a in chemicalName.ChemicalName.to_list():
    chemList.append(a)

#Setting Shape attribute to note Nodes that are the original chemicals to change shape for clarity
chemsShape={chem: 'diamond' for chem in chemList}
nx.set_node_attributes(G, chemsShape,'shape')

#Setting Size attribute to note Nodes that are the original chemicals to change size for clarity
chemSize={chem: 30 for chem in chemList}
nx.set_node_attributes(G, chemSize,'size')
##



##Figure out the proper layout for the graph
net = Network(width="1000px", height="700px", bgcolor='#222222', font_color='white',  select_menu=True, filter_menu=True)
pos = nx.spring_layout(G,scale=5000)
net.from_nx(G)

#Turn off physics
net.toggle_physics(False)
for node in net.get_nodes():
    net.get_node(node)['x']=pos[node][0]
    net.get_node(node)['y']=-pos[node][1] #the minus is needed here to respect networkx y-axis convention 
    net.get_node(node)['physics']=False
    net.get_node(node)['label']=str(node) #set the node label so that it can be displayed

#write and save a file of the graph
net.save_graph(name = "Pyvis_Graph.html")
print("Pyvis_Graph.html saved. That's the visualization and thats the end of the program/workflow.\nHave fun researching!")