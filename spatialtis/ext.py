import shutil
from ast import literal_eval
from pathlib import Path
from typing import Union

import pandas as pd
from anndata import AnnData

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc


@doc
class prepare_svca(AnalysisBase):
    """Prepare data for SVCA analysis

    Spatial Variance Components Analysis: `SVCA <https://github.com/damienArnol/svca>`_

    The input format is separated folder for each ROI with `expressions.txt` and `positions.txt`.

    Args:
        data: {adata}
        out_dir: The directory to store the data
        entry: The name of new folder to store the data
        **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        out_dir: Union[Path, str],
        entry: str = "svca_data",
        **kwargs
    ):
        super().__init__(data, task_name="prepare_svca", **kwargs)
        groups = data.obs.groupby(self.exp_obs)

        p = Path(out_dir) / entry

        if p.exists():
            shutil.rmtree(p)

        p.mkdir(exist_ok=True)

        for n, g in groups:
            folder_name = "_".join([str(i) for i in n])
            folder_path = p / folder_name
            folder_path.mkdir()

            expression_path = folder_path / "expressions.txt"
            position_path = folder_path / "positions.txt"

            # export to expressions.txt
            pd.DataFrame(
                data.X[[int(i) for i in g.index]], columns=data.var[self.marker_key]
            ).to_csv(expression_path, sep="\t", index=False)

            # export to positions.txt
            cents = []
            for c in g[self.centroid_key]:
                c = literal_eval(c)
                cents.append([c[0], c[1]])
            pd.DataFrame(cents).to_csv(
                position_path, sep="\t", index=False, header=False
            )
