import pandas as pd
import os
from bokeh.plotting import figure, show, from_networkx
from bokeh.models import ColumnDataSource, LinearColorMapper
from sklearn.neighbors import radius_neighbors_graph
import networkx as nx
from bokeh.models import Circle

def load_data(panel_type, datafolder):
    files = [f for f in os.listdir(datafolder) if f.endswith(f'{panel_type}.csv')]
    data_frames = [pd.read_csv(f'{datafolder}/{file}').assign(source_file=file) for file in files]
    data = pd.concat(data_frames)
    return data

option = 'IF1'  # Option can be 'IF1', 'IF2', or 'IF3'
datafolder = 'if_data_short'
data = load_data(option, datafolder)[:400000] #data reduced because of bad performance
print('Data loaded')

def create_phenotype_graph(data, phenotype, radius):
    phenotype_data = data[data['phenotype'].str.contains(f'{phenotype}\+')]
    coordinates = phenotype_data[['nucleus.x', 'nucleus.y']].values
    graph = radius_neighbors_graph(coordinates, radius=radius, mode='connectivity', include_self=False)

    # Create a NetworkX graph from the sparse matrix
    G = nx.from_scipy_sparse_array(graph)

    # Create a Bokeh plot
    plot = figure(title="Interactive Graph Visualization", x_range=(-1, 1), y_range=(-1, 1),
                  tools="", toolbar_location=None)

    graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))
    plot.renderers.append(graph_renderer)

    show(plot)

create_phenotype_graph(data, 'CD20', 30)
