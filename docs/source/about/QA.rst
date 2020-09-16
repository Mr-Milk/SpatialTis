Q&A
====

Basic
------

What is spatialtis for?
#######################

SpatialTis is a spatial analysis package for single cell multiplexed data.

How parallel processing works?
#################################

The parallel happens at ROI level. If a small dataset is used, there might find the processing time will take longer than non-parallel mode, because there are some overhead steps (serialization and deserialization so that it could distribute the data to different processors). This features only support Linux and MacOS, on Windows there will be no effects.

Only following functions in spatialtis have parallel processing support, pass argument `mp=True` to enable.

    - :class:`spatialtis.read_ROIs`
    - :func:`spatialtis.spatial.spatial_distribution`
    - :func:`spatialtis.spatial.hotspot`

These are implemented in rust, it will automatically run in parallel.

    - :meth:`spatialtis.Neighbors.find_neighbors`
    - :func:`spatialtis.spatial.neighborhood_analysis`
    - :func:`spatialtis.spatial.spatial_enrichment_analysis`

MacOS Issues
-------------

Do you want the application "Python.app" to accept the incoming network connection?
#####################################################################################

.. image:: ../src/mac_connections_issue.png
    :width: 50%
    :align: center

If there are lots of these windows pop up on Mac, it's cause by *Ray*.
If you find it annoying, the simplest solution is to turn off your firework (with safety risk)
or add it into the firewall white list.
Another solution from this `stackoverflow answer <https://stackoverflow.com/a/59186900>`_ might also be helpful.
