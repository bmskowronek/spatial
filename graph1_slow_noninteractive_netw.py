import pandas as pd
from sklearn.neighbors import radius_neighbors_graph
import os
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def load_data(panel_type, datafolder):
    files = [f for f in os.listdir(datafolder) if f.endswith(f'{panel_type}.csv')]
    data_frames = [pd.read_csv(f'{datafolder}/{file}').assign(source_file=file) for file in files]
    data = pd.concat(data_frames)
    return data

option = 'IF1'  # Option can be 'IF1', 'IF2', or 'IF3'
datafolder = 'if_data_short'
data = load_data(option, datafolder)[:60000]
#print(data[:60])
print('Data loaded')


def create_phenotype_graph(data, phenotype, radius):
    '''
    output: modified phenotype_data dataframe which for each cell with requested phenotype
            has an entry in the column "connections" which is a tuple with indexes of
            cells that it is connected to

    '''
    # Filter data for the specified phenotype
    phenotype_data = data[data['phenotype'].str.contains(f'{phenotype}\+')]
    # Extract coordinates
    coordinates = phenotype_data[['nucleus.x', 'nucleus.y']].values
    # Create the radius neighbors graph
    graph = radius_neighbors_graph(coordinates, radius=radius, mode='connectivity', include_self=False)
    
    rows, cols = graph.nonzero()
    graph_tuples = list(zip(rows, cols))
    # Create a defaultdict to hold the groups
    grouped = defaultdict(list)
    
    # Populate the defaultdict with grouped values
    for k, v in graph_tuples:
        grouped[k].append(v)
    
    # Determine the maximum index to ensure coverage of all possible indices
    max_index = len(coordinates)
    
    # Prepare the result list with handling for skipped indices
    result = []
    for i in range(max_index):
        if i in grouped:
            result.append(tuple(grouped[i]))
        else:
            result.append((None,))  # Add a tuple with None for skipped indices
    phenotype_data = phenotype_data.assign(connections=result)
    phenotype_data = phenotype_data.reset_index(drop=True)
    return phenotype_data
    

graph_data = create_phenotype_graph(data, 'CD20', 30)

G = nx.Graph()

# Add nodes
for index, row in graph_data.iterrows():
    print(index)
    G.add_node(index, color=row['color'])

# Add edges
for index, row in graph_data.iterrows():
    if row['connections'] != (None,):
        connections = row['connections']  # Convert string tuple to actual tuple
        for connection in connections:
            if connection in data.index:  # Check if the connection index exists in the data
                G.add_edge(index, connection)
                print(index)

# Extract colors for corresponding nodes
colors = [G.nodes[node]['color'] for node in G.nodes()]


# Assuming 'cell type' and 'color' are columns in your graph_data DataFrame
unique_cell_types = graph_data[['cell type', 'color']].drop_duplicates()
legend_patches = [mpatches.Patch(color=row['color'], label=row['cell type']) for index, row in unique_cell_types.iterrows()]
plt.figure(figsize=(16, 16))
nx.draw(G, node_color=colors, with_labels=True, node_size=500, font_size=4, font_color='black')
print('plot created')
plt.title('Graph Visualization of Cells')

# Add legend to the plot
plt.legend(handles=legend_patches, title="Cell Types", loc='best')
plt.show()
