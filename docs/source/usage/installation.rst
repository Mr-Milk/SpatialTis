Installation
============

SpatialTis requires **Python >= 3.6**, it's recommended that you install it in a new environment.

.. note::
    Make sure you use SpatialTis with >= 0.4.0, I moved most of the implementation to Rust with huge
    performance improve. Most of the parallelization don't rely on ray anymore. And I dropped the interactive
    visualization support.

pypi
----
Install the basic of spatialtis::

    pip install spatialtis

For the full features::

    pip install spatialtis[all]

    # In some terminal environment you may try
    pip install 'spatialtis[all]'

