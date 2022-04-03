FROM pytorchlightning/pytorch_lightning
WORKDIR /work
RUN pip install --no-cache-dir scanpy spatialtis[all]
RUN pip install --no-cache-dir jupyter jupyterlab jupyter_http_over_ws dask-labextension
RUN jupyter serverextension enable --py jupyter_http_over_ws dask_labextension && ipcluster nbextension disable
WORKDIR /analysis
CMD ["jupyter", "lab", "--notebook-dir=/analysis", "--port=8888", "--ip=0.0.0.0", "--allow-root", "--no-browser"]
