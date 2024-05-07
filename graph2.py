import pandas as pd
from sklearn.neighbors import radius_neighbors_graph
import os
from collections import defaultdict



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

import plotly.graph_objects as go

def plot_interactive_graph(data):
    # Create a scatter plot of the cell positions
    trace_points = go.Scatter(
        x=data['nucleus.x'],
        y=data['nucleus.y'],
        mode='markers',
        marker=dict(size=5, color='blue'),
        name='Cells',
        text=data['phenotype'],
        hoverinfo='text'
    )

    # Prepare to plot the connections
    edge_x = []
    edge_y = []
    
    # Iterate over the DataFrame to plot lines between connected nodes
    for idx, row in data.iterrows():
        if row['connections'][0] != (None,):  # Ensure there are connections
            for conn_idx in row['connections']:
                if conn_idx is not None:
                    # Start point
                    edge_x.append(row['nucleus.x'])
                    edge_y.append(row['nucleus.y'])
                    # End point
                    edge_x.append(data.loc[conn_idx, 'nucleus.x'])
                    edge_y.append(data.loc[conn_idx, 'nucleus.y'])
                    # Add None to x, y lists to prevent continuous lines
                    edge_x.append(None)
                    edge_y.append(None)

    trace_lines = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode='lines',
        line=dict(width=1, color='grey'),
        hoverinfo='none'
    )

    # Define layout
    layout = go.Layout(
        title='Cell Connectivity Graph',
        title_x=0.5,
        xaxis=dict(title='Nucleus X'),
        yaxis=dict(title='Nucleus Y'),
        showlegend=False
    )

    # Create figure and add traces
    fig = go.Figure(data=[trace_lines, trace_points], layout=layout)
    fig.show()

# Call the function with the graph data
plot_interactive_graph(graph_data)



