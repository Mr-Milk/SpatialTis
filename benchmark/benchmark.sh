# python fake_data.py
#python plot.py

docker run -it --rm -d --name spatialtis_1c8g --mount type=bind,source="$PWD",target=/usr/src/app/benchmark --cpus 1 --memory 8g time_spatialtis bash
docker exec -it spatialtis_1c8g /bin/sh -c "cd benchmark; python time_spatialtis.py 1"
docker stop spatialtis_1c8g


docker run -it --rm -d --name spatialtis_4c16g --mount type=bind,source="$PWD",target=/usr/src/app/benchmark --cpus 4 --memory 16g time_spatialtis bash
docker exec -it spatialtis_4c16g /bin/sh -c "cd benchmark; python time_spatialtis.py 4"
docker stop spatialtis_4c16g


docker run -it --rm -d --name spatialtis_8c32g --mount type=bind,source="$PWD",target=/usr/src/app/benchmark --cpus 8 --memory 32g time_spatialtis bash
docker exec -it spatialtis_8c32g /bin/sh -c "cd benchmark; python time_spatialtis.py 8"
docker stop spatialtis_8c32g

docker run -it --rm -d --name spatialtis_12c64g --mount type=bind,source="$PWD",target=/usr/src/app/benchmark --cpus 12 --memory 64g time_spatialtis bash
docker exec -it spatialtis_12c64g /bin/sh -c "cd benchmark; python time_spatialtis.py 12"
docker stop spatialtis_12c64g

docker run -it --rm -d --name squidpy_12c64g --mount type=bind,source="$PWD",target=/usr/src/app/benchmark --cpus 12 --memory 64g time_squidpy bash
docker exec -it squidpy_12c64g /bin/sh -c "cd benchmark; python time_squidpy.py 12"
docker stop squidpy_12c64g


docker run -it --rm -d --name giotto_1c8g --mount type=bind,source="$PWD",target=/home/your_giotto/benchmark --cpus 1 --memory 8g delron01/giotto bash
docker exec -it giotto_1c8g /bin/sh -c "cd benchmark; Rscript time_giotto.r 1"
docker stop giotto_1c8g


docker run -it --rm -d --name giotto_4c16g --mount type=bind,source="$PWD",target=/home/your_giotto/benchmark --cpus 4 --memory 16g delron01/giotto bash
docker exec -it giotto_4c16g /bin/sh -c "cd benchmark; Rscript time_giotto.r 4"
docker stop giotto_4c16g


docker run -it --rm -d --name giotto_8c64g --mount type=bind,source="$PWD",target=/home/your_giotto/benchmark --cpus 8 --memory 64g delron01/giotto bash
docker exec -it giotto_8c48g /bin/sh -c "cd benchmark; Rscript time_giotto.r 8"
docker stop giotto_8c48g

rm rprof.log