Save Results
=======================

Save analysis results
++++++++++++++++++++++

There are two methods to save the analysis results of SptialTis.

Just directly access :code:`.result` attribute::

    result = st.cell_components(data).result

Or you can read it from anndata if you know the key::

    result = st.get_result(data, 'your_key')


Save visualization results
++++++++++++++++++++++++++++

If you set the :code:`CONFIG.AUTO_SAVE`, the image can be saved automatically::

    CONFIG.AUTO_SAVE = True # create a default saved path
    CONFIG.AUTO_SAVE = "your_saved_folder" # provide your own saved path

Or you can manually save it::

    sp.cell_components(data).save("saved_path/img.png")

The default will ensure a publication quality save (dpi>=300).
To gain more access to the plot, you can get the plot instance::

    sp.cell_components(data).fig # For matplotlib
    sp.cell_components(data).plot # For bokeh & pyecharts

.. note::
    Starting from 0.3.0, SpatialTis provides static visualization based on matplotlib for every analysis. The interactive
    visualizations are saved to `.html` by default. If you want to export them into static images, the best way is to
    use screen capture.

All the interactive visualization are rendered in HTML5 Canvas/SVG. To save these plots into static
images, extra dependencies need to be installed. There are three visualization libraries used in spatialtis:

+--------------+----------------------+-----------------------+-------------------+
|              | **bokeh**            | **pyecharts**         | **matplotlib**    |
+--------------+----------------------+-----------------------+-------------------+
| Renderer     | Canvas,              | Canvas,               | WX, GTK,          |
|              | SVG                  | SVG (Not perfect)     | Cairo ...         |
+--------------+----------------------+-----------------------+-------------------+
| Naive export | `.html`, `.svg`      | `.html`, `.png`       | `.png`, `.jpeg`,  |
|              |                      |                       | `.eps`, `.pdf` ...|
+--------------+----------------------+-----------------------+-------------------+


To export `.svg` in pyecharts, pass parameter `renderer='svg'` into plotting function to switch to svg backend,
the default is 'canvas'; However, the svg is not perfect, your exported image might has layout issue.
This can only be solved if echarts.js support `svg` further in the future.


For bokeh, you need to install `geckodriver <https://github.com/mozilla/geckodriver/releases>`_ with
`firefox <https://www.mozilla.org/firefox/new/>`_
or `chromedriver <https://chromedriver.chromium.org/downloads>`_ with `chromium <https://download-chromium.appspot.com/>`_,
check `selenium-python <https://selenium-python.readthedocs.io/installation.html#drivers>`_, download selenium afterwards;
For pyecharts, `phantom.js` is recommended.

Installation of selenium
--------------------------

- `selenium-python <https://selenium-python.readthedocs.io/installation.html#drivers>`_

::

    # use conda
    conda install -c conda-forge selenium
    # use pip
    pip install selenium

Installation of phantomjs
--------------------------
::

    # use conda
    conda install -c conda-forge phantomjs
    # use npm
    npm install phantomjs

Installation of firefox and geckodriver
----------------------------------------------------

- `geckodriver <https://github.com/mozilla/geckodriver/releases>`_
- `firefox <https://www.mozilla.org/firefox/new/>`_

::

    conda install -c conda-forge firefox geckodriver

Or you can download them from the official release site, remember to add geckodriver to your path if
install manually.

Installation of chromium and chromedriver
------------------------------------------

- `chromedriver <https://chromedriver.chromium.org/downloads>`_
- `chromium <https://download-chromium.appspot.com/>`_

If you have chrome on you system, you only need to install chromedriver, remember to add chromedriver to your path.
