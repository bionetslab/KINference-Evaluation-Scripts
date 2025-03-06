import pandas as pd
from arboreto.algo import grnboost2
from distributed import LocalCluster, Client

# Load the data
if __name__ == '__main__':
    x0 = pd.read_csv('../data/Bouhaddou2023/intensities_Mock10h.tsv', sep = '\t', index_col = 0)
    x0.dropna(inplace=True)
    x0 = x0.T
    x1 = pd.read_csv('../data/Bouhaddou2023/intensities_VIC10h.tsv', sep = '\t', index_col = 0)
    x1.dropna(inplace=True)
    x1 = x1.T

    local_cluster = LocalCluster(n_workers=30,
                                 threads_per_worker=1,
                                 memory_limit=8e9)
    custom_client = Client(local_cluster)

    # Build custom "tf" (in this case all the kinase substrates) list
    serine_threonine_kinases = pd.read_csv('../data/kinase_data/serine_threonine_kinases/kinase_name_mappings.tsv', sep = '\t', index_col = 0).loc[:, 'ACC#'].unique()
    tyrosine_kinases = pd.read_csv('../data/tyrosine_kinases/kinase_name_mappings.tsv', sep = '\t', index_col = 0).loc[:, 'ACC#'].unique()

    kinases = list(serine_threonine_kinases) + list(tyrosine_kinases)
    tf_names = [kin for kin in x0.columns if str.split(kin, '_')[0] in kinases]

    # x0 network
    print(f'Building GRNBoost2 network for x0')
    network = grnboost2(expression_data = x0, tf_names = tf_names, client_or_address = custom_client)
    network.to_csv('../results/GRNBoost2/x0_network.tsv', sep = '\t', index = False, header = False)

    # x1 network
    tf_names = [kin for kin in x1.columns if str.split(kin, '_')[0] in kinases]
    print(f'Building GRNBoost2 network for x1')
    network = grnboost2(expression_data = x1, tf_names = tf_names, client_or_address = custom_client)
    network.to_csv('../results/GRNBoost2/x1_network.tsv', sep = '\t', index = False, header = False)

    custom_client.close()
    local_cluster.close()
