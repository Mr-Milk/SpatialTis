Rscript time_giotto.r
python time_spatialtis.py
python plot.py


docker run -it --rm \
 --mount type=bind,source=/Users/milk/Documents/project/neighborhood_performance_compare,target=/usr/src/app \
 --cpus 1 \
 time_spatialtis bash

docker run -it --rm \
 --mount type=bind,source=/Users/milk/Documents/project/neighborhood_performance_compare,target=/usr/src/app \
 --cpus 1 \
 time_giotto bash