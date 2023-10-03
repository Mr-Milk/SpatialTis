FROM python:3.12
WORKDIR /usr/src/app
RUN pip install --no-cache-dir squidpy memory_profiler