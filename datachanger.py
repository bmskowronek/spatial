import streamlit as st
import pandas as pd
import os
import bisect


def rgb_to_hex(rgb):
    '''
    expected input is like: (190, 255, 255) etc
    expected output is a hex color: '#beffff' or whatever
    '''
    return '#%02x%02x%02x' % rgb


IF1_cell_mapping = {"other": rgb_to_hex((190, 190, 190)),
                    "CD15+Tumor": rgb_to_hex((73, 176, 248)),
                    "CD15-Tumor": rgb_to_hex((138, 79, 45)),
                    "Tcell": rgb_to_hex((235, 74, 148)),
                    "Bcell": rgb_to_hex((204, 49, 31)),
                    "BnTcell": rgb_to_hex((236, 95, 42)),
                    "Neutrophil": rgb_to_hex((0, 40, 245)),
                    "Macrophage": rgb_to_hex((97, 209, 62)),
                    "DC": rgb_to_hex((49, 113, 30))}

IF2_cell_mapping = {"Epithelial cells": rgb_to_hex((138, 79, 45)),
                    "CD8+_Tcells TOTAL": rgb_to_hex((235, 74, 148)),
                    "Ki67+_CD8+_Tcells": rgb_to_hex((236, 95, 42)),
                    "Ki67+_CK+": rgb_to_hex((73, 176, 211)),
                    "CK+_Ki67+_PDL1+": rgb_to_hex((49, 113, 30)),
                    "CK+_PDL1+": rgb_to_hex((97, 209, 62)),
                    "PD1+": rgb_to_hex((204, 49, 31)),
                    "other": rgb_to_hex((190, 190, 190)),
                    'CD8+ T cell':rgb_to_hex((180, 102, 25)),
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
                    "CD8+/CD4+ T cell": rgb_to_hex((236, 155,2)),
                    "FOXP3+ Treg":rgb_to_hex((255, 127, 127)),
                    "Tumor": rgb_to_hex((111, 34, 4)),
                    "other": rgb_to_hex((190, 190, 190)),
                    "NK cell":rgb_to_hex((141, 95, 176)),
                    "not_defined":rgb_to_hex((77, 70, 82))}
cell_dict = {
    "IF1": IF1_cell_mapping,
    "IF2": IF2_cell_mapping,
    "IF3": IF3_cell_mapping,}

def pheno_binarify(phenotype):
     '''
     in: string 'CD8-CK+GB+Ki67-PD1+PDL1+'
     out: list sorted alphabetically ['CD8-', 'CK+', 'GB+', 'Ki67-', 'PD1+', 'PDL1+']
            and then turned into a string of binary characters
     the genes are sorted to maintain some sort of logic behind their ordering,
     so that it stays consistent

        last line turns ['CD8-', 'CK+', 'GB+', 'Ki67-', 'PD1+', 'PDL1+']
        into '011011'
        elements with - become 0s, + become 1s, it can be consistent due to alphabetical sorting
     '''
     found_genes = []
     i=0
     j=0
     while i<len(phenotype):
         if phenotype[i] not in '+-':
             i += 1
             j += 1
         else:
             gene = phenotype[(i-j):i+1]
             bisect.insort(found_genes, gene) #faster sorting by inserting like this... maybe?
             j = 0
             i += 1
     return ''.join(['0' if '-' in item else '1' for item in found_genes])




# Path to the folder containing the CSV files
folder_path = 'if_data'
folder_path2 = 'if_data_short'

# Get the list of CSV files and group them based on the end of the file name
if_files = {}
for filename in os.listdir(folder_path):
    if filename.endswith(".csv"):
        if filename.endswith("IF1.csv"):
            if_files.setdefault("IF1", []).append(filename)
        elif filename.endswith("IF2.csv"):
            if_files.setdefault("IF2", []).append(filename)
        elif filename.endswith("IF3.csv"):
            if_files.setdefault("IF3", []).append(filename)

def pheno_to_celltype(phenotype):
    binary_phenotype = pheno_binarify(phenotype)
    cell_type = phenotype_to_celltype_dict[binary_phenotype]
    return cell_type

def if_coloring(celltype, panelname):
    if panelname == 'IF1':
        return IF1_cell_mapping[celltype]
    elif panelname== 'IF2':
        return IF2_cell_mapping[celltype]
    elif panelname == 'IF3':
        return IF3_cell_mapping[celltype]
    else:
        print('unknown panel name chosen, try IF1, IF2, IF3')

chosen_panel = 'IF3'
colors = pd.read_csv(f"{chosen_panel}_phen_to_cell_mapping.csv")
phenotype_to_celltype_dict = dict(zip([pheno_binarify(col) for col in colors['phenotype']], colors['celltype']))
selected_files=if_files[chosen_panel]
for filename in selected_files:
        data = pd.read_csv(f'if_data//{filename}')
        phenotypes = data['phenotype']
        data['cell type'] = data['phenotype'].apply(pheno_to_celltype)
        data['color'] = data['cell type'].apply(if_coloring, panelname= chosen_panel)
        data_small = data[['cell.ID', 'nucleus.x', 'nucleus.y', 'phenotype', 'cell type', 'color']]
        data_small.to_csv(f'{folder_path2}//{str(filename)}', sep=',', index=False, encoding='utf-8')
