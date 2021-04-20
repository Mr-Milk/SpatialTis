FROM giotto
WORKDIR /usr/src/app
COPY Giotto-1.0.3.tar.gz ./
RUN R -e "remotes::install_github('RubD/Giotto');"