import pandas as pd
import numpy as np
import argparse
import pcst_fast

# defining parser
parser = argparse.ArgumentParser()
parser.add_argument('--net_file', help='Path to file with connected baseline KIN')
parser.add_argument('--node_names_file', help='Path to file with nodenames')
parser.add_argument('--output_path', help='Path to output folder')

# parse arguments
args = parser.parse_args()

# PCST
edges_df = pd.read_csv(args.net_file, sep = '\t')
nodes_df = pd.read_csv(args.node_names_file, sep = '\t')

e = edges_df[['Source_idx', 'Target_idx']]
node_prizes = np.absolute(nodes_df['f']).to_numpy()

# defining pcst algorithm parameters
root = -1
num_clusters = 1
pruning = 'strong'
verbosity_level = 0
diameter = 6

# setting edges costs uniformly to 1
edge_costs = np.ones(e.shape[0])

vertices, edges = pcst_fast.pcst_fast(e, node_prizes, edge_costs, root, num_clusters, pruning, verbosity_level)

steiner_edges = edges_df.iloc[edges, :]
steiner_edges.to_csv(args.output_path, sep = '\t', index = False)