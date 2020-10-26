Installation
============

SpatialTis requires **Python >= 3.6**, it's recommended that you install it in a new environment.

pypi
----
Install the basic of spatialtis::

    pip install spatialtis

For the full features::

    pip install spatialtis[all]

    # In some terminal environment you may try
    pip install 'spatialtis[all]'

For windows user
-----------------

Using `WSL2 on Windows 10 <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_
might relief you from misery installation process.

If run on Windows, you may find installation problem with following dependencies, please refer to their documentation for detail instructions.
You can try install them via conda channels, or compile it from source on your machine:

- `python-igraph <https://igraph.org/python/>`_
- `leidenalg <https://leidenalg.readthedocs.io/en/stable/install.html>`_
- `ray <https://docs.ray.io/en/latest/installation.html>`_: You need to install the `Visual C++ Runtime <https://aka.ms/vs/16/release/vc_redist.x64.exe>`_

Complete installation
----------------------
If you want to save the visualization of SpatialTis to static images,
extra dependencies are needed, please check `Export to static image <image_export.html>`_.

If you work in **Zepplin** or **Jupyter Lab**, please check `Working in notebook <notebook.html>`_.