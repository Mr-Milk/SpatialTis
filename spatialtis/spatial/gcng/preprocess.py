from collections import Counter

import numpy as np
import pandas as pd


def neighbors_pairs(labels, neighbors):
    p1, p2 = [], []
    for l, ns in zip(labels, neighbors):
        for n in ns:
            if n > l:
                p1.append(l)
                p2.append(n)
    return np.array([p1, p2], dtype=np.int32)


def overlap_genes(data_genes, lr_genes):
    overlap_set = set([i for i in data_genes]).intersection(set([i for i in lr_genes]))
    return list(overlap_set)


def balanced_pairs(pairs):
    """Generate the same number of negative pairs as positive pairs"""
    pairs.iloc[:, 2] = pairs.iloc[:, 2].astype(dtype=float)
    existed_pairs = [tuple(i) for i in pairs.iloc[:, [0, 1]].to_numpy()]
    relationship = pairs.iloc[:, 2].values
    count = Counter(relationship)
    add_amount = count[1] - count.get(0, 0)

    if add_amount > 0:
        pos_pairs = pairs[pairs.iloc[:, 2] == 1]
        p1 = pos_pairs.iloc[:, 0].values
        p2 = pos_pairs.iloc[:, 1].values

        neg_count = 0
        neg_pairs = []
        while neg_count < count[1]:
            ix1, ix2 = np.random.randint(0, len(p1), 2)
            p = (p1[ix1], p2[ix2])
            if p not in existed_pairs:
                neg_pairs.append((p1[ix1], p2[ix2], 0.0))
                neg_count += 1

        return pd.concat([pairs, pd.DataFrame(neg_pairs, columns=pairs.columns)]).reset_index(drop=True)
    else:
        return pairs


def train_test_split(pairs, partition=0.9):
    """Split the train and test dataset, ligand or receptor appear in train will not show in test
    The negative pairs will be randomly generated, the number of negative pairs will be equal to positive pairs

    """
    pairs = pairs.copy()
    pairs.iloc[:, 2] = pairs.iloc[:, 2].astype(dtype=float)
    ix = np.arange(0, len(pairs))
    cut = int(len(ix) * partition)
    train = pairs.iloc[:cut, :]
    test = pairs.iloc[cut:, :]
    train_pool = np.unique(train.iloc[:, [0, 1]].to_numpy())
    clean_test = []
    for i, row in test.iterrows():
        p1, p2, _ = row
        if (p1 not in train_pool) & (p2 not in train_pool):
            clean_test.append(row)
    test = pd.DataFrame(clean_test)
    test = balanced_pairs(test)
    train = balanced_pairs(train)
    return train, test


def graph_data_loader(pairs, exp, markers_mapper, npairs, device="cpu", batch_size=32, shuffle=True):
    try:
        import torch
        from torch_geometric.data import Data
        from torch_geometric.loader import DataLoader
    except ImportError as e:
        raise ImportError("Please install pytorch and torch_geometric")

    dataset = []
    edge_index = torch.tensor(npairs, dtype=torch.long)
    for _, row in pairs.iterrows():
        ix1 = int(markers_mapper[row[0]])
        ix2 = int(markers_mapper[row[1]])

        x = torch.tensor([exp[ix1], exp[ix2]], dtype=torch.float).T
        t = torch.tensor([row[2]], dtype=torch.float)
        data = Data(x=x, edge_index=edge_index, y=t)
        data.to(device)
        dataset.append(data)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=True)


def predict_data_loader(pairs, exp, markers_mapper, npairs, device="cpu", batch_size=32):
    try:
        import torch
        from torch_geometric.data import Data
        from torch_geometric.loader import DataLoader
    except ImportError as e:
        raise ImportError("Please install pytorch and torch_geometric")

    dataset = []
    edge_index = torch.tensor(npairs, dtype=torch.long)
    for _, row in pairs.iterrows():
        ix1 = int(markers_mapper[row[0]])
        ix2 = int(markers_mapper[row[1]])

        x = torch.tensor([exp[ix1], exp[ix2]], dtype=torch.float).T
        data = Data(x=x, edge_index=edge_index)
        data.to(device)
        dataset.append(data)
    return DataLoader(dataset, batch_size=batch_size)
