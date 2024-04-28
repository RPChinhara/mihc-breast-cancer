import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import os
import math
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")


def mapping(df):
    mapping = {
        'CD3+CD4':'Helper T cell',
        'CD8':'Killer T cell',
        'CD20':'B cell',
        'CD3':'T cell',
        'CD68':'Macrophage',
        'CD3+CD8':'Killer T cell',
        'CD4':'Helper T cell',
        'FoxP3':'Regulatory T cell',
        'CD3+FoxP3':'Regulatory T cell',
        'CD4+FoxP3': 'Regulatory T cell',
        'CD8+FoxP3': 'Regulatory T cell',
        'CD3+CD4+FoxP3': 'Regulatory T cell',
        'CD3+CD8+FoxP3': 'Regulatory T cell',
        'CD3+CD4+CD8':'CD4+CD8'
    }   

    df['Class'] = df['Class'].map(lambda x: mapping.get(x, f'{x}'))

    return df

def data_cleaning(df):
    df = df[["Class", "Parent", "Centroid X µm", "Centroid Y µm"]]
    df['Class'] = df['Class'].fillna(df['Parent'])
    df = df.drop(['Parent'], axis=1)

    df['Class'] = df['Class'].apply(lambda x: x.replace('Target: ', ''))
    df['Class'] = df['Class'].apply(lambda x: x.replace(': ', '+'))

    new_rows = {
        'Class': [], 
        'Centroid X µm': [],
        'Centroid Y µm': []
    }

    for idx, row in df.iterrows():
        if 'CD20' in list(row['Class'].split('+')) or 'CD68' in list(row['Class'].split('+')):
            for clas_ in list(row['Class'].split('+')):
                new_rows['Class'].append(clas_)

            for _ in range(len(list(row['Class'].split('+')))):
                new_rows['Centroid X µm'].append(row['Centroid X µm'])
                new_rows['Centroid Y µm'].append(row['Centroid Y µm']) 
            
            df = df.drop(idx)  # type: ignore

    new_rows = pd.DataFrame(new_rows)
    df = pd.concat([df, new_rows], ignore_index=True) 

    df = mapping(df)

    new_rows = {
        'Class': [], 
        'Centroid X µm': [],
        'Centroid Y µm': []
    }

    for idx, row in df.iterrows():
        if 'CD4' in list(row['Class'].split('+')) and 'CD8' in list(row['Class'].split('+')):
            for clas_ in list(row['Class'].split('+')):
                new_rows['Class'].append(clas_)

            for _ in range(len(list(row['Class'].split('+')))):
                new_rows['Centroid X µm'].append(row['Centroid X µm'])
                new_rows['Centroid Y µm'].append(row['Centroid Y µm']) 
            
            df = df.drop(idx)  # type: ignore

    new_rows = pd.DataFrame(new_rows)
    df = pd.concat([df, new_rows], ignore_index=True)

    df = mapping(df)

    return df

def euclidean_distance(p1, p2):
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def network_viz(df, fname):
    df = data_cleaning(df)

    G = nx.Graph()
    for idx, row in df.iterrows():
        G.add_node(idx, pos=(row['Centroid X µm'], - row['Centroid Y µm']), type=row['Class'])

    for (node1, row1), (node2, row2) in combinations(df.iterrows(), 2):
        pos1 = (row1['Centroid X µm'], row1['Centroid Y µm'])
        pos2 = (row2['Centroid X µm'], row2['Centroid Y µm'])
        if euclidean_distance(pos1, pos2) < 35:
            G.add_edge(node1, node2)

    nx.write_gml(G, fname)

    # tumor_nodes = [node for node, data in G.nodes(data=True) if data.get('type') == 'Tumor']
    # tumor_graph = G.subgraph(tumor_nodes)

    # pos = nx.spring_layout(tumor_graph)

    # plt.figure(figsize=(6, 6))
    # nx.draw(tumor_graph, pos, node_size=10, node_color='steelblue')
    # plt.savefig(fname=fname, pad_inches=0)


if __name__=='__main__':
    folder_path = 'rem_patients_data/'

    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    for csv_file in csv_files:
        if csv_file == 'Detections_OP_P5.csv':
            patient_id = csv_file[-6:-4]
            file_path = os.path.join(folder_path, csv_file)
            df = pd.read_csv(file_path)
            network_viz(df, f'rem_patients_gml/P{patient_id}.gml')
            break
