import pandas as pd
import os
from collections import defaultdict
import plotly.graph_objects as go
import numpy as np
from sklearn.neighbors import radius_neighbors_graph


def load_data(panel_type, datafolder):
    files = [f for f in os.listdir(datafolder) if f.endswith(f'{panel_type}.csv')]
    data_frames = [pd.read_csv(f'{datafolder}/{file}').assign(source_file=file) for file in files]
    data = pd.concat(data_frames)
    return data


option = 'IF1'  # Option can be 'IF1', 'IF2', or 'IF3'
datafolder = 'if_data_short'
data = load_data(option, datafolder)
print('Data loaded')


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
    print(len(phenotype_data))
    usable_data = phenotype_data[phenotype_data.connections != (None,)]
    print(len(usable_data))
    n = 2
    usable_data = usable_data[usable_data['connections'].apply(lambda x: len(x) > n)] # only cells with more than n neighbors
    print(len(usable_data))

    fig = go.Figure()

    grouped = usable_data.groupby('source_file')

    for name, group in sorted(grouped):
        fig.add_trace(
            go.Scattergl(
                mode='markers',
                x=group['nucleus.x'],
                y=group['nucleus.y'],
                marker=dict(
                    color=group['color'],
                    size=10  # Increased marker size
                ),
                name=name,
                hovertemplate=(
                        "Cell Type: %{customdata[0]}<br>" +
                        "Cell ID: %{customdata[1]}<br>" +
                        "Source File: %{customdata[2]}<br>" +
                        "Phenotype:%{customdata[2]}<br>" +
                        "<extra></extra>"
                ),
                customdata=np.stack((group['cell type'], group['cell.ID'], group['source_file'],group['phenotype']), axis=-1),
                showlegend=True
            )
        )

        # Add lines between connected cells


    fig.update_layout(
        title=f'Interactive scatter plot of cell data - phenotypes {phenotype}+, neighbor radius {radius}, showing only cells with more than {n} neighbors',
        width=1600,
        height=900,
        template='plotly_dark'
    )

    fig.show()
    fig.write_html(f"plot{option}_neighbors.html")
    print('plot loaded')

create_phenotype_graph(data, 'CD20', 30)

