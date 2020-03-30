Export to static image
=======================

.. note::
    If you need to run on cluster or non-display environment, it's recommended that you save the data first and visualize it
    in your own computer with screen. For pyecharts' visualization, you should export it as html, open in browser, adjust it
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
| Other export | `.png`, `.pdf`,      | `.pdf`, `.svg`,       |                   |
|              | `.eps`               | `.eps`                |                   |
+--------------+----------------------+-----------------------+-------------------+

All static export will ensure the image quality reach the normal publication requirement (dpi>=300), if you need further
configuration, use `return_plot=True` to get the plot instance.


For bokeh >= 2.0.0, you need to install `geckodriver` with `firefox` or `chromedriver` with `chromium browser`,
check `selenium-python <https://selenium-python.readthedocs.io/installation.html#drivers>`_, download selenium afterwards;
For pyecharts, `phantom.js` is recommended::

    # install via conda
    conda install selenium phantomjs

    # install selenium via pip
    pip install selenium

    # install phantomjs via npm/cnpm
    npm install phantomjs
    # if you are in China, pull from Alibaba's mirror
    cnpm install phantomjs

To export `.svg` in pyecharts, use `renderer='svg'` to switch to svg backend, the default is 'canvas'; However, the svg
is not perfect, your exported image might has layout issue. This can only be solved if echarts.js support `svg`
further in the future.
