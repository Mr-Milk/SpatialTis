# python fake_data.py
#python plot.py

docker run -it --rm -d \
  --name spatialtis_1c8g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 1 \
  --memory 8g \
  time_spatialtis bash
docker exec -it spatialtis_1c8g /bin/sh -c "cd benchmark; python time_spatialtis.py 1"
docker stop spatialtis_1c8g


docker run -it --rm -d \
  --name spatialtis_4c16g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 4 \
  --memory 16g \
  time_spatialtis bash
docker exec -it spatialtis_4c16g /bin/sh -c "cd benchmark; python time_spatialtis.py 4"
docker stop spatialtis_4c16g


docker run -it --rm -d \
  --name spatialtis_8c32g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 8 \
  --memory 32g \
  time_spatialtis bash
docker exec -it spatialtis_8c32g /bin/sh -c "cd benchmark; python time_spatialtis.py 8"
docker stop spatialtis_8c32g


docker run -it --rm -d \
  --name giotto_1c8g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 1 \
  --memory 8g \
  time_giotto bash
docker exec -it giotto_1c8g /bin/sh -c "cd benchmark; Rscript time_giotto.r 1"
docker stop giotto_1c8g


docker run -it --rm -d \
  --name giotto_4c16g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 4 \
  --memory 16g \
  time_giotto bash
docker exec -it giotto_4c16g /bin/sh -c "cd benchmark; Rscript time_giotto.r 4"
docker stop giotto_4c16g


docker run -it --rm -d \
  --name giotto_8c32g \
  --mount type=bind,source="$PWD",target=/usr/src/app/benchmark \
  --cpus 8 \
  --memory 32g \
  time_giotto bash
docker exec -it giotto_8c32g /bin/sh -c "cd benchmark; Rscript time_giotto.r 8"
docker stop giotto_8c32g

rm rprof.log