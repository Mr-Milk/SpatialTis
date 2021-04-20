#/usr/bin sh

cd ..
echo "$PWD"
docker run -it --rm -d \
  --name spatialtis_4c16g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 4 \
  --memory 16g \
  time_spatialtis bash
docker exec -it spatialtis_4c16g /bin/sh -c "cd benchmark/benchmark; python spatialtis_full_analysis.py"
docker stop spatialtis_4c16g
