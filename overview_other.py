import pandas as pd
import os
import plotly.graph_objects as go
import numpy as np

# Function to convert RGB to hex
def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

# Color mappings for different panels
IF1_cell_mapping = {
    "other": rgb_to_hex((190, 190, 190)),
    "CD15+Tumor": rgb_to_hex((73, 176, 248)),
    "CD15-Tumor": rgb_to_hex((138, 79, 45)),
    "Tcell": rgb_to_hex((235, 74, 148)),
    "Bcell": rgb_to_hex((204, 49, 31)),
    "BnTcell": rgb_to_hex((236, 95, 42)),
    "Neutrophil": rgb_to_hex((0, 40, 245)),
    "Macrophage": rgb_to_hex((97, 209, 62)),
    "DC": rgb_to_hex((49, 113, 30))
}
IF2_cell_mapping = {"Epithelial cells": rgb_to_hex((138, 79, 45)),
                    "CD8+_Tcells TOTAL": rgb_to_hex((235, 74, 148)),
                    "Ki67+_CD8+_Tcells": rgb_to_hex((236, 95, 42)),
                    "Ki67+_CK+": rgb_to_hex((73, 176, 211)),
                    "CK+_Ki67+_PDL1+": rgb_to_hex((49, 113, 30)),
                    "CK+_PDL1+": rgb_to_hex((97, 209, 62)),
                    "PD1+": rgb_to_hex((204, 49, 31)),
                    "other": rgb_to_hex((190, 190, 190)),
                    "CD8+ T cell": rgb_to_hex((235, 74, 148)),
                    "CD8+/CD4+ T cell": rgb_to_hex((236, 155,2)),
                    "CD4+ T cell": rgb_to_hex((204, 49, 31)),
                    "FOXP3+ Treg":rgb_to_hex((255, 127, 127)),
                    "NK cell":rgb_to_hex((141, 95, 176)),
                    "T cell": rgb_to_hex((236, 95, 42)),
                    "not_defined":rgb_to_hex((77, 70, 82)),
                    "Tumor": rgb_to_hex((73, 176, 248))}

IF3_cell_mapping = {"CD8+ T cell": rgb_to_hex((235, 74, 148)),
                    "T cell": rgb_to_hex((236, 95, 42)),
                    "CD4+ T cell": rgb_to_hex((204, 49, 31)),
                    "CD8+/CD4+ T cell": rgb_to_hex((255, 5, 113)),
                    "FOXP3+ Treg":rgb_to_hex((255, 127, 127)),
                    "Tumor": rgb_to_hex((111, 34, 4)),
                    "other": rgb_to_hex((190, 190, 190)),
                    "NK cell":rgb_to_hex((141, 95, 176)),
                    "not_defined":rgb_to_hex((77, 70, 82))}

merged_cell_mapping = {**IF2_cell_mapping, **IF3_cell_mapping, **IF1_cell_mapping}

print(merged_cell_mapping)
# Function to load data
datafolder = 'if_data_short'

def load_data(panel_type, datafolder):
    files = [f for f in os.listdir(datafolder) if f.endswith(f'{panel_type}.csv')]
    data_frames = [pd.read_csv(f'{datafolder}/{file}').assign(source_file=file) for file in files]
    data = pd.concat(data_frames)
    return data

option = 'IF3'  # Option can be 'IF1', 'IF2', or 'IF3'

data = load_data(option, datafolder)
print('Data loaded')

# Count the occurrences of each cell type
cell_type_counts = data['cell type'].value_counts()

# Calculate the percentage of each cell type
total_counts = cell_type_counts.sum()
cell_type_percentages = (cell_type_counts / total_counts * 100).round(2)



# Create a histogram of the cell type counts with colors and percentage labels
fig = go.Figure(data=[go.Bar(
    x=cell_type_counts.index,
    y=cell_type_counts.values,
    text=cell_type_percentages.apply(lambda x: f'{x}%'),
    textposition='outside',
    marker_color=[merged_cell_mapping.get(cell, '#000000') for cell in cell_type_counts.index]  # Use color mapping
)])

fig.update_layout(
    title=f'Histogram of Cell Types - {option}',
    xaxis_title='Cell Type',
    yaxis_title='Count',
    autosize=False,
    width=1000,
    height=600
)

fig.show()