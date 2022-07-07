Installation
============

SpatialTis requires **Python >= 3.8**.

.. |Version Support| image:: https://img.shields.io/pypi/pyversions/spatialtis?style=flat-square

PYPI
----
Install the basic of spatialtis

.. code-block:: bash

    pip install spatialtis

For the full features

.. code-block:: bash

    pip install 'spatialtis[all]'


Docker
-------
The quickest way to run is to use a docker image.

.. code-block:: bash

    docker pull mr-milk/spatialtis

To run a jupyter notebook from the docker image and mount your data folder to it:

.. code-block:: bash

    cd your/data/
    docker run -it --rm -p 8888:8888 -v "${PWD}:/analysis" spatialtis
