import logging
from itertools import product
from pathlib import Path
from typing import Optional, Union, List, Tuple

import pandas as pd
from anndata import AnnData
from spatialtis.abc import AnalysisBase
from spatialtis.utils import read_neighbors, log_print, doc

from .preprocess import overlap_genes, train_test_split, neighbors_pairs, \
    graph_data_loader, predict_data_loader

MODEL_SAVE_KEY = "SpatialTis-GCNG-Model-State"


@doc
class GCNG(AnalysisBase):
    """A pytorch reimplementation of GCNG

    Use to identify directional gene-gene interactions. The trained model will be automatically save to anndata.

    .. note::
        To perform this analysis, you need `pytorch <https://pytorch.org/get-started/locally/#start-locally>`_,
        `pytorch-geometry <https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html>`_ and
        `pytorch-lightning <https://www.pytorchlightning.ai/>`_ installed. If you have GPU, make sure you install
        pytorch with GPU support, it would be way more faster than CPU.

    Args:
        data: {adata}
        known_pairs: The input data for training, should be a dataframe with three columns,
            ligand, receptor, relationship; 0 means not interact, 1 means interact.
        predict_pairs: The pairs that you interested
        train_partition: The ratio to split the dataset for training
        gpus: Number of gpu to use, can be auto-detected
        max_epochs: Number of epoch
        lr: Learning rate
        batch_size: The batch size
        random_seed: The random seed
        load_model: To load a pretrained model from anndata
        **kwargs: {analysis_kwargs}

    """

    def __init__(
            self,
            data: AnnData,
            known_pairs: Optional[pd.DataFrame] = None,
            predict_pairs: Optional[List[Tuple]] = None,
            train_partition: float = 0.9,
            gpus: Optional[int] = None,
            max_epochs: int = 10,
            lr: float = 1e-4,
            batch_size: int = 32,
            random_seed: int = 42,
            load_model: bool = False,
            **kwargs,
    ):
        try:
            import torch
            import torch.nn.functional as F
            from torch.nn import Flatten, Linear
            from torch_geometric.nn import GCNConv, global_max_pool
            import pytorch_lightning as pl
            from pytorch_lightning.core.lightning import LightningModule
            logging.getLogger("lightning").setLevel(0)
        except ImportError:
            raise ImportError("To run GCNG, please install pytorch, pytorch-lightning, "
                              "torch-geometric, torch_sparse and torch_scatter.")
        super().__init__(data, display_name="GCNG", **kwargs)
        device = "cpu"
        if gpus is None:
            cuda_count = torch.cuda.device_count()
            gpus = cuda_count
            if gpus > 0:
                device = "cuda"

        # To make pytorch a optional deps
        # We could only init the model from within
        class GCNGModel(LightningModule):
            def __init__(self, node_size, output_features, lr=lr):
                super().__init__()
                self.conv1 = GCNConv(2, 32)
                self.conv2 = GCNConv(32, 32)
                self.dense1 = Linear(output_features * node_size * 32, 512)
                self.dense2 = Linear(512, output_features)
                self.flatten = Flatten()

                self.lr = lr
                self.correct = 0
                self.test_data_len = 0
                self.acc = 0
                self.pred = []

            def forward(self, x, edge_index):
                # x, edge_index = data.x, data.edge_index
                x = self.conv1(x, edge_index)
                x = F.elu(x)
                x = self.conv2(x, edge_index)
                x = F.elu(x)
                x = torch.flatten(x)
                x = self.dense1(x)
                x = F.elu(x)
                x = self.dense2(x)

                return torch.sigmoid(x)

            def configure_optimizers(self):
                return torch.optim.Adam(self.parameters(), lr=self.lr)

            def training_step(self, train_data, batch_idx):
                x, edge_index, batch = train_data.x, train_data.edge_index, train_data.batch
                x = self(x, edge_index)
                loss_in = x.flatten()
                loss_out = train_data.y
                loss = F.binary_cross_entropy(loss_in, loss_out)
                return loss

            def test_step(self, test_data, batch_idx):
                x, edge_index, batch = test_data.x, test_data.edge_index, test_data.batch
                x = self(x, edge_index)
                pred = x.detach().cpu().numpy().flatten().round()
                truth_y = test_data.y.cpu().numpy()
                self.correct += (pred == truth_y).sum()
                self.test_data_len += len(test_data.y)
                self.acc = self.correct / self.test_data_len
                return self.acc

            def predict_step(self, predict_data, batch_idx):
                x, edge_index, batch = predict_data.x, predict_data.edge_index, predict_data.batch
                x = self(x, edge_index)
                self.pred += x.detach().cpu().numpy().flatten().round().tolist()
                return self.pred

        # init model and trainer first
        gc = GCNGModel(data.n_obs, batch_size, lr=lr)
        pl.seed_everything(random_seed, workers=True)
        trainer = pl.Trainer(gpus=gpus,
                             max_epochs=max_epochs,
                             deterministic=True,
                             progress_bar_refresh_rate=0,
                             weights_summary=None)
        # create neighbors pairs
        npairs = neighbors_pairs(data.obs[self.cell_id_key],
                                 read_neighbors(data.obs, self.neighbors_key))
        # get exp info and create markers mapper
        # markers' name will all be lowercase
        exp = data.X.T
        markers = data.var.markers
        markers = pd.Series([i.lower() for i in markers], index=markers.index)
        markers_mapper = dict(zip(markers.tolist(), range(len(markers))))

        if load_model:  # load pre-trained model
            try:
                state = self.data.uns[MODEL_SAVE_KEY]
            except KeyError:
                raise ValueError("Pre-trained model not found, please retrain the model")
            gc.load_state_dict(state)
        else:  # train the model
            # find overlap genes
            lr_genes = known_pairs.iloc[:, [0, 1]].to_numpy().flatten()
            overlap_sets = overlap_genes(self.markers, lr_genes)
            filtered_pairs = known_pairs[known_pairs.iloc[:, 0].isin(overlap_sets) &
                                         known_pairs.iloc[:, 1].isin(overlap_sets)].iloc[:, [0, 1, 2]]
            if len(filtered_pairs) == 0:
                raise ValueError("The gene in `known_pairs` has no overlap with genes in data")
            train, test = train_test_split(filtered_pairs, train_partition)
            train_loader = graph_data_loader(train, exp, markers_mapper, npairs, device, batch_size, shuffle=True)
            test_loader = graph_data_loader(test, exp, markers_mapper, npairs, device, batch_size, shuffle=False)
            trainer.fit(gc, train_loader)
            trainer.test(dataloaders=test_loader, verbose=False)
            log_print(f"Model accuracy {gc.acc}")
            self.data.uns[MODEL_SAVE_KEY] = gc.state_dict()  # save model
            self.model = gc  # allow user to access model
            if predict_pairs is not None:
                pairs_all = set(["*".join(p) for p in product(markers, repeat=2)])
            else:
                pairs_all = predict_pairs
            pairs_db = set(["*".join(p) for p in filtered_pairs.iloc[:, [0, 1]].to_numpy()])
            predict_pairs = [p.split("*") for p in pairs_all.difference(pairs_db)]
        # the model output is dynamically adjust according to batch size
        # the predict step should be able to iter through all pairs
        if predict_pairs is None:
            raise ValueError("To rerun the model, you must specific the `predict_pairs`"
                             "and tell spatialtis the ligand-receptor pairs you want to predict.")
        predict_size = len(predict_pairs)
        append_amount = batch_size - predict_size % batch_size
        predict_pairs += predict_pairs[:append_amount]
        predict = pd.DataFrame(predict_pairs)
        predict_loader = predict_data_loader(predict, exp, markers_mapper, npairs, device, batch_size)
        # init the model and train
        trainer.predict(dataloaders=predict_loader)
        predict['relationship'] = gc.pred
        predict.columns = ['Gene1', 'Gene2', 'relationship']
        self.result = predict.iloc[:predict_size, :].copy()