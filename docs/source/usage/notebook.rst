Working in notebook
--------------------

If you work with other notebook environments that's not jupyter notebook, you can config it by ::

    from spatialtis import CONFIG
    CONFIG.WORKING_ENV = 'zepplin'
    # ["jupyter_notebook", "jupyter_lab", "nteract", "zepplin", None]

If working in Zepplin, please install `bkzep`::

    pip install bkzep

If working in Jupyter Lab, please install the following dependencies::

    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter labextension install @bokeh/jupyter_bokeh

If working in headless mode, you should just save the plot.