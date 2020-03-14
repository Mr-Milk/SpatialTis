Q&A
====

Basic
------

What is spatialtis for?
#######################

SpatialTis

Introduction to these technologies
###################################

If you want to know about these technologies,

How parallel processing works?
#################################

The parallel happens at ROI level. If a small dataset is used, there might find the processing time will take longer than non-parallel mode, because there are some overhead steps (serialization and deserialization so that it could distribute the data to different processors). This features only support Linux and MacOS, on Windows there will be no effects.

Only following functions in spatialtis have parallel processing support, pass argument `mp=True` to enable.

    - `spatialtis.read_ROIs <api_index/preprocessing.rst#spatialtis.read_ROIs>`_
    - spatialtis.Neighbors.find_neighbors()
    - spatialtis.spatial.neighborhood_analysis()
    - spatialtis.spatial.spatial_enrichment_analysis()
    - spatialtis.spatial.hotspot()


Some plotting function not showing images?
###########################################

Some of spatialtis visualization use bokeh, which is rendered in canvas, if you work with other notebook environments not jupyter, bokeh supports zepplin, you can config it by ::

    from spatialtis import CONFIG

    CONFIG.WORKING_ENV = 'zepplin'

If working in Jupyter Lab/Hub, please install the following dependencies::

    jupyter labextension install @bokeh/jupyter_bokeh

If working with Zeppelin, please install `bkzep` (Only PyPI)::

    pip install bkzep

If working in headless mode, you should just save the plot.

MacOS Issues
-------------

Do you want the application "Python.app" to accept the incoming network connection?
#####################################################################################

.. image:: ../src/mac_connections_issue.png
    :width: 50%
    :align: center

If there are lots of these windows pop up on Mac, it's cause by *Ray*. If you find it annoying, the simplest solution is to turn off your firework (with safety risk) or add it into the firewall white list. Another solution from this `stackoverflow answer <https://stackoverflow.com/a/59186900>`_ might also be helpful.
