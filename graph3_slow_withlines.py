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
data = load_data(option, datafolder)[:300000] #data reduced because of bad performance
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

    fig = go.Figure()

    # Add scatter plot for cells
    grouped = phenotype_data.groupby('source_file')
    for name, group in grouped:
        fig.add_trace(
            go.Scattergl(
                mode='markers',
                x=group['nucleus.x'],
                y=group['nucleus.y'],
                marker=dict(
                    color=group['color'],
                    size=10
                ),
                name=name,
                hovertemplate=(
                    "Cell Type: %{customdata[0]}<br>" +
                    "Cell ID: %{customdata[1]}<br>" +
                    "Source File: %{customdata[2]}<br>" +
                    "Phenotype: %{customdata[3]}<br>"
                    "<extra></extra>"
                ),
                customdata=np.stack((group['cell type'], group['cell.ID'], group['source_file'], group['phenotype']), axis=-1)
            )
        )

    # Add lines for connections
    for index, row in phenotype_data.iterrows():
        if row['connections'] != (None,):
            for connected_index in row['connections']:
                if connected_index is not None:
                    connected_cell = phenotype_data.iloc[connected_index]
                    fig.add_trace(
                        go.Scattergl(
                            x=[row['nucleus.x'], connected_cell['nucleus.x']],
                            y=[row['nucleus.y'], connected_cell['nucleus.y']],
                            mode='lines',
                            line=dict(color='gray', width=1),
                            showlegend=False
                        )
                    )

    fig.update_layout(
        title=f'Interactive scatter plot of cell data - phenotypes {phenotype}+',
        width=1600,
        height=900,
        template='plotly_dark'
    )

    fig.write_html(f"plot{option}graph_test.html")
    print('plot loaded')


    
create_phenotype_graph(data, 'CD20', 30)
datatime2 = time.time()
print('\n time is: '+ str(datatime2-start))

