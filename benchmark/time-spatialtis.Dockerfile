FROM python:3
WORKDIR /usr/src/app
RUN pip install --no-cache-dir numpy && pip install --no-cache-dir spatialtis memory_profiler
