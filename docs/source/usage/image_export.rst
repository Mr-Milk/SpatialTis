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
|              | SVG                  | SVG (No perfect)      | Cairo ...         |
+--------------+----------------------+-----------------------+-------------------+
| Naive export | `.html`, `.svg`      | `.html`, `.png`       | `.png`, `.jpeg`,  |
|              |                      |                       | `.eps`, `.pdf` ...|
+--------------+----------------------+-----------------------+-------------------+


All static export will ensure the image quality reach the normal publication requirement (dpi>=300), if you need further
configuration, use `return_plot=True` to get the plot instance.


For bokeh >= 2.0.0, you need to install `geckodriver <https://github.com/mozilla/geckodriver/releases>`_ with
`firefox <https://www.mozilla.org/firefox/new/>`_
or `chromedriver <https://chromedriver.chromium.org/downloads>`_ with `chromium <https://download-chromium.appspot.com/>`_,
check `selenium-python <https://selenium-python.readthedocs.io/installation.html#drivers>`_, download selenium afterwards;
For pyecharts, `phantom.js` is recommended::

    # Option 1: install via conda
    conda install -c conda-forge selenium phantomjs
    conda install -c conda-forge geckodriver firefox

    # Option 2: install selenium via pip and phantomjs via npm
    pip install selenium
    npm install phantomjs


To export `.svg` in pyecharts, pass parameter `renderer='svg'` into plotting function to switch to svg backend,
the default is 'canvas'; However, the svg
is not perfect, your exported image might has layout issue. This can only be solved if echarts.js support `svg`
further in the future.
