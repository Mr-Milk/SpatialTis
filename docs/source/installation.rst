Installation
============

spatialTis requires **Python 3.6+**, it has full features support for MacOS and Linux.

For windows users:

Parallelism is not available for windows and will not be supported.
If you use Windows 10 version 1903 or higher, you can try WSL to use full features of spatialTis.
Docker is another solution.


Anaconda
---------
We highly recommend you install via anaconda, if you haven't install anaconda, please check [link to conda]

First by create a new conda virtual environment::

    conda create -n tissue-analysis # the name of the env is tissue-analysis

And then activate the environment::

    conda activate tissue-analysis

Now we can install it::

    conda install -c bioconda spatialtis


pypi
----
::

    pip install spatialtis


Docker
------
We provide a mini docker image, it contains some useful python packages besides spatialtis.

Start jupyter notebook from docker
