#!/usr/bin/env/ python
#Here we want to make a script which translates our graph to the same format as the LastGraph which is the output format of Velvet

import networkx as nx
import json
from networkx.readwrite import json_graph
import http_server


def Convert():

	#Frist read in the graph
	G = nx.read_gpickle("test.gpickle")

	######################
	# Convert to json    #
	######################

	#this d3 example uses the name attribute for the mouse-hover value, so add a name to each node
	for n in G:
		G.node[n]['name'] = len(G.node[n]['Base'])
		#G.node[n]['name'] = G.node[n]['Base']

	#write json formatted data
	d = json_graph.node_link_data(G)

	#write json
	json.dump(d,open('force/force.json','w'))
	print('Wrote node-link JSON data to force/force.json')
	
	#open URL in running web browser
	http_server.load_url('force/index2.html')
	print('Or copy all files in force/ to webserver and load force/force.html')

if __name__ == '__main__':
	Convert()
