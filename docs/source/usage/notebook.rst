Working in notebook
--------------------

If you work with other notebook environments that's not jupyter notebook, interactive visualization might not render properly.
You need manually set::

    # For pyecharts
    from pyecharts.globals import CurrentConfig, NotebookType
    CurrentConfig.NOTEBOOK_TYPE = NotebookType.JUPYTER_LAB # Using JupyterLab
    CurrentConfig.NOTEBOOK_TYPE = NotebookType.NTERACT # Using Nteract
    CurrentConfig.NOTEBOOK_TYPE = NotebookType.ZEPPELIN # Using Zeppelin


If working in Zepplin, please install `bkzep` to get `bokeh` running::

    pip install bkzep

If working in Jupyter Lab, please install the following dependencies::

    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter labextension install @bokeh/jupyter_bokeh

If working in headless mode, you should just save the plot.