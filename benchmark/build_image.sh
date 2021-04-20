mkdir SpatialTis
cp -r ../spatialtis SpatialTis
cp ../setup.py SpatialTis
cp ../README.md SpatialTis
docker build -t time_spatialtis -f time-spatialtis.Dockerfile .
rm -rf ./SpatialTis

# docker build -f time-giotto.Dockerfile .