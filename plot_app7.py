import streamlit as st
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


# Streamlit user interface
st.title('Interactive Scatter Plot of Cell Data')

# Dropdown to select the panel type
option = st.selectbox(
    'Choose a panel type:',
    ('IF1', 'IF2', 'IF3')
)

# Button to trigger plot creation
if st.button('Create Plot'):
    data = load_data(option, datafolder)
    print('data loaded')
    # Create a Plotly scatter plot
    x=np.array(data['nucleus.x'])
    y=np.array(data['nucleus.y'])
    colors = np.array(data['color'])
    #info for dot labels
    cell_types = np.array(data['cell type'])
    cell_ids = np.array(data['cell.ID'])
    source_files = np.array(data['source_file'])
    # Build figure
    fig = go.Figure()

    # Add scatter trace
    fig.add_trace(
        go.Scattergl(
            mode='markers',
            x=x,
            y=y,
            marker=dict(
                color=colors,
                size=5
            ),
            hovertemplate=(
                    "Cell Type: %{customdata[0]}<br>" +
                    "Cell ID: %{customdata[1]}<br>" +
                    "Source File: %{customdata[2]}<br>" +
                    "<extra></extra>"
            ),
            customdata=np.stack((cell_types, cell_ids, source_files), axis=-1),
            showlegend=False
        )
    )
    fig.update_layout(
        width=1600,
        height=900,
        template='plotly_dark'
    )
    # Display the plot
    st.plotly_chart(fig, use_container_width=True)
    print('plot loaded')
    fig.write_html(f"plot{option}.html")
    st.write("Plot creation complete!")
