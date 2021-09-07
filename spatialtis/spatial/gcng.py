from typing import List, Tuple

from anndata import AnnData

from spatialtis.abc import AnalysisBase


class GCNG(AnalysisBase):
    def __init__(self, data: AnnData,
                 known_pairs: List[Tuple[str, str]],
                 ):
        super().__init__(data)

        # create positive pairs / negative pairs
