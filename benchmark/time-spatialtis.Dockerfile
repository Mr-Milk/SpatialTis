FROM rayproject/ray
WORKDIR /usr/src/app
COPY ./SpatialTis/ .
RUN pip install -e . && pip install memory_profiler
