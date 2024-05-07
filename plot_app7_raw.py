import pandas as pd
import os
import plotly.graph_objects as go
import numpy as np

# Function to load data
datafolder = 'if_data_short'

def load_data(panel_type, datafolder):
    # List comprehension to find all files ending with the specified panel_type
    files = [f for f in os.listdir(datafolder) if f.endswith(f'{panel_type}.csv')]

    # Create a list of DataFrames, each with a new column 'source_file' indicating the source file name
    data_frames = [pd.read_csv(f'{datafolder}/{file}').assign(source_file=file) for file in files]

    # Concatenate all DataFrames into a single DataFrame
    data = pd.concat(data_frames)

    return data

option = 'IF1'

data = load_data(option, datafolder)
print('data loaded')

# Create a Plotly scatter plot
fig = go.Figure()

# Group data by source file
grouped = data.groupby('source_file')

# Iterate over each group to create a scatter trace
for name, group in sorted(grouped):
    fig.add_trace(
        go.Scattergl(
            mode='markers',
            x=group['nucleus.x'],
            y=group['nucleus.y'],
            marker=dict(
                color=group['color'],
                size=5
            ),
            name=name,  # Use the file name as the trace name
            hovertemplate=(
                "Cell Type: %{customdata[0]}<br>" +
                "Cell ID: %{customdata[1]}<br>" +
                "Source File: %{customdata[2]}<br>" +
                "<extra></extra>"
            ),
            customdata=np.stack((group['cell type'], group['cell.ID'], group['source_file']), axis=-1),
            showlegend=True
        )
    )

fig.update_layout(
    title='Interactive scatter plot of cell data',
    width=1600,
    height=900,
    template='plotly_dark'
)

fig.write_html(f"plot{option}.html")
print('plot loaded')