import networkx as nx
from statistics import mean
import os
import sys
import pandas as pd

def get_mixing_score(G):
    number_of_immune_cells = sum(1 for node in G.nodes() if G.nodes[node]['type'] in ['Helper T cell', 'Killer T cell', 'T cell', 'B cell', 'Macrophage'])

    if number_of_immune_cells < 250:
        return -1

    immune_to_tumor = sum(1 for u, v in G.edges() if (G.nodes[u]['type'] == 'Tumor' and G.nodes[v]['type'] in ['Helper T cell', 'Killer T cell', 'T cell', 'B cell', 'Macrophage']) or (G.nodes[v]['type'] == 'Tumor' and G.nodes[u]['type'] in ['Helper T cell', 'Killer T cell', 'T cell', 'B cell', 'Macrophage']))

    immune_to_immune = sum(1 for u, v in G.edges() if (G.nodes[u]['type'] in ['Helper T cell', 'Killer T cell', 'T cell', 'B cell', 'Macrophage'] and G.nodes[v]['type'] in ['Helper T cell', 'Killer T cell', 'T cell', 'B cell', 'Macrophage']))

    return immune_to_tumor / immune_to_immune

def categorize_mixing_score(G):
    mixing_score = get_mixing_score(G)
    if mixing_score == -1:
        return (mixing_score, 'Cold')
    elif mixing_score < 0.22:
        return (mixing_score, 'Compartmentalized')
    else:
        return (mixing_score, 'Mixed')

def get_stromal_clustering(G):
    stroma_subgraph = G.subgraph([node for node in G.nodes() if G.nodes[node]['type'] == 'Other'])
    avg_stroma_cc = nx.average_clustering(stroma_subgraph, nodes=[node for node in G.nodes() if G.nodes[node]['type'] == 'Other'])
    return avg_stroma_cc

def get_stromal_barrier(G, node_type):
    b_cells = [node for node, data in G.nodes(data=True) if data.get('type') == node_type]
    tumor_nodes = [node for node, data in G.nodes(data=True) if data.get('type') == 'Tumor']

    shortest_paths = {}
    closest_tumor = None
    closest_path = 1000000

    for b_cell in b_cells:
        shortest_paths[b_cell] = {}
        
        for tumor_node in tumor_nodes:
            try:
                length = nx.shortest_path_length(G, source=b_cell, target=tumor_node)
            except:
                continue

            if length < closest_path and length > 1:
                    closest_path = length
                    closest_tumor = tumor_node
        
        shortest_paths[b_cell] = {closest_tumor: closest_path}

    num_other_nodes = {}
    for b_cell, paths in shortest_paths.items():
        num_other_nodes[b_cell] = {}
        for tumor_node, length in paths.items():
            path = nx.shortest_path(G, source=b_cell, target=tumor_node)
            num_other = sum(1 for node in path if G.nodes[node]['type'] == 'Other')
            num_other_nodes[b_cell][tumor_node] = num_other

    avg_num_other_nodes = sum(sum(num_other_nodes[b_cell].values()) for b_cell in num_other_nodes) / len(num_other_nodes)

    return avg_num_other_nodes

if __name__=='__main__':
    file_path = sys.argv[1]

    print(file_path)

    G = nx.read_gml(file_path)

    mixing_score = categorize_mixing_score(G)

    stromal_cc = get_stromal_clustering(G)

    stromal_barrier_data = []
    for node_type in ['B cell', 'T cell', 'Macrophage', 'Killer T cell', 'Helper T cell', 'Regulatory T cell']:
        barrier = get_stromal_barrier(G, node_type=node_type)
        stromal_barrier_data.append(barrier)

    data = pd.DataFrame({
        'Mixing Score': [mixing_score[0]],
        'Phenotype': [mixing_score[1]],
        'Stromal Clustering': [stromal_cc],
        'Stromal Barrier - B cell': [stromal_barrier_data[0]],
        'Stromal Barrier - T cell': [stromal_barrier_data[1]],
        'Stromal Barrier - Macrophage': [stromal_barrier_data[2]],
        'Stromal Barrier - Killer T cell': [stromal_barrier_data[3]],
        'Stromal Barrier - Helper T cell': [stromal_barrier_data[4]],
        'Stromal Barrier - Regulatory T cell': [stromal_barrier_data[5]]
    })

    data.to_csv(f'P{file_path[-6:-4]}.csv')