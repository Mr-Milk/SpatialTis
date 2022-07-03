FROM python:3.10
RUN pip install --no-cache-dir scanpy spatialtis==0.5.0rc0 lightgbm scikit-image leidenalg \
    jupyter jupyterlab jupyter_http_over_ws ipywidgets
WORKDIR /analysis
EXPOSE 8888
CMD ["jupyter", "lab", "--notebook-dir=/analysis", "--port=8888", "--ip=0.0.0.0", "--allow-root", "--no-browser"]
