import pandas as pd
import os
from collections import defaultdict
import plotly.graph_objects as go
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
import time

start = time.time()
def load_data(panel_type, datafolder):
    files = [f for f in os.listdir(datafolder) if f.endswith(f'{panel_type}.csv')]
    data_frames = [pd.read_csv(f'{datafolder}/{file}').assign(source_file=file) for file in files]
    data = pd.concat(data_frames)
    return data


option = 'IF1'  # Option can be 'IF1', 'IF2', or 'IF3'
datafolder = 'if_data_short'
data = load_data(option, datafolder)
print(len(data))
print('Data loaded')
datatime = time.time()
print('\n time is: '+ str(datatime-start))

def create_phenotype_graph(data, phenotype, radius):
    phenotype_data = data[data['phenotype'].str.contains(f'{phenotype}\+')]
    coordinates = phenotype_data[['nucleus.x', 'nucleus.y']].values
    graph = radius_neighbors_graph(coordinates, radius=radius, mode='connectivity', include_self=False)

    rows, cols = graph.nonzero()
    graph_tuples = list(zip(rows, cols))
    grouped = defaultdict(list)

    for k, v in graph_tuples:
        grouped[k].append(v)

    max_index = len(coordinates)
    result = []
    for i in range(max_index):
        if i in grouped:
            result.append(tuple(grouped[i]))
        else:
            result.append((None,))
    phenotype_data = phenotype_data.assign(connections=result)
    phenotype_data = phenotype_data.reset_index(drop=True)
    return phenotype_data

dax = create_phenotype_graph(data, 'CD20', 30)

shorter =  dax[['nucleus.x', 'nucleus.y', 'connections']].copy()



class KDNode:
    def __init__(self, x, y, index, connections, left=None, right=None):
        self.x = x  # The x coordinate
        self.y = y  # The y coordinate
        self.index = index  # The index of the point in the dataset
        self.connections = connections  # List of connected indexes
        self.left = left  # Left child
        self.right = right  # Right child

class KDTree:
    def __init__(self):
        self.root = None

    def insert(self, x, y, index, connections, depth=0):
        def insert_rec(node, x, y, index, connections, depth):
            if node is None:
                return KDNode(x, y, index, connections)

            # Calculate current dimension (x=0, y=1)
            cd = depth % 2

            if (x if cd == 0 else y) < (node.x if cd == 0 else node.y):
                node.left = insert_rec(node.left, x, y, index, connections, depth + 1)
            else:
                node.right = insert_rec(node.right, x, y, index, connections, depth + 1)

            return node

        self.root = insert_rec(self.root, x, y, index, connections, depth)

    def print_tree(self, node, depth=0):
        if node is not None:
            self.print_tree(node.left, depth + 1)
            print(f"{'  ' * depth}Index: {node.index}, X: {node.x}, Y: {node.y}, Connections: {node.connections}")
            self.print_tree(node.right, depth + 1)

# Example usage
kd_tree = KDTree()


# Iterate over DataFrame rows as tuples
for row in shorter.itertuples(index=True, name=None):
    index = row[0]
    x = row[1]
    y = row[2]
    connections = row[3]
    kd_tree.insert(x, y, index, connections)

# Optionally, print the KD-Tree
#kd_tree.print_tree(kd_tree.root)























