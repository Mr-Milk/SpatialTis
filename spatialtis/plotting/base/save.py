# Save function for bokeh and pyecharts
from pathlib import Path

from bokeh.io import export_png, export_svgs, output_file, save
from pyecharts.render import make_snapshot
from snapshot_phantomjs import snapshot


def save_bokeh(plot, path):
    save_path = Path(path)
    file_ext = save_path.suffix[1:]
    save_path = str(save_path)

    plot.background_fill_color = None
    plot.border_fill_color = None

    if file_ext not in ["html", "svg", "png"]:
        raise NotImplementedError("Current supported formats are: svg, html and png")
    elif file_ext == "html":
        output_file(save_path)
        save(plot)
    elif file_ext == "svg":
        plot.output_backend = "svg"
        export_svgs(plot, filename=save_path)
    else:
        export_png(plot, filename=save_path)


def save_pyecharts(charts, path, delay=2, pixel_ratio=5):
    save_path = Path(path)
    file_ext = save_path.suffix[1:]
    charts.renderer = "canvas"

    if file_ext not in ["html", "png"]:
        raise ValueError("Current supported formats are: png and html")
    elif file_ext == "html":
        charts.render(save_path)
    else:
        make_snapshot(
            snapshot,
            charts.render(),
            str(save_path),
            delay=delay,
            pixel_ratio=pixel_ratio,
        )
