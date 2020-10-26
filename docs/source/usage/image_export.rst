Export to static image
=======================

.. note::
    If you need to run on cluster or non-display environment, it's recommended that you save the data first and visualize it
    in your own computer with screen. So that you can have better authetic control towards the visuzliation.
    For pyecharts' visualization, you should export it as html, open in browser, adjust it
    till you are satisfied and then use the save bottom to save it in `.png`. If you run in notebook, bokeh and pyecharts can
    be rendered inline.

Some of the visualizations in spatialtis are interactive, it's rendered in HTML5 Canvas/SVG. To save these plots into static
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


All static export will ensure the image quality reach the normal publication requirement (dpi>=300), if you need further
configuration, use `return_plot=True` to get the plot instance.

To export `.svg` in pyecharts, pass parameter `renderer='svg'` into plotting function to switch to svg backend,
the default is 'canvas'; However, the svg
is not perfect, your exported image might has layout issue. This can only be solved if echarts.js support `svg`
further in the future.


For bokeh >= 2.0.0, you need to install `geckodriver <https://github.com/mozilla/geckodriver/releases>`_ with
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
