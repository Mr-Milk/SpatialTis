# Save function for bokeh and pyecharts
from pathlib import Path
from bokeh.io import export_svgs, output_file, save
from cairosvg import svg2pdf, svg2png, svg2ps

from pyecharts.render import make_snapshot
from snapshot_phantomjs import snapshot


def save_bokeh(plot, path, dpi=300):
    save_path = Path(path)
    file_ext = save_path.suffix[1:]

    if file_ext not in ['html', 'svg', 'png', 'pdf', 'ps']:
        raise ValueError("File extension must be specified, supported formats are: svg, html and (png, pdf, ps)")
    elif file_ext is 'html':
        output_file(save_path)
        save(plot)
    elif file_ext is 'svg':
        plot.output_backend = "svg"
        export_svgs(plot, filename=save_path)
    else:
        # for other formats, it will be saved as .svg and then convert to other format
        plot.output_backend = "svg"
        tmp_svg = str(save_path) + ".svg"
        export_svgs(plot, filename=tmp_svg)

        if file_ext is 'png':
            svg2png(url=tmp_svg, write_to=save_path, dpi=dpi)
        elif file_ext is 'pdf':
            svg2pdf(url=tmp_svg, write_to=save_path, dpi=dpi)
        elif file_ext is 'ps':
            svg2ps(url=tmp_svg, write_to=save_path, dpi=dpi)


def save_pyecharts(charts, path, delay=2, pixel_ratio=5):
    save_path = Path(path)
    file_ext = save_path.suffix[1:]
    charts.renderer = 'canvas'

    if file_ext not in ['html', 'svg', 'png', 'pdf', 'eps']:
        raise ValueError("File extension must be specified, supported formats are: png, html and (svg, pdf, eps)")
    elif file_ext is 'html':
        charts.render(save_path)
    elif file_ext is 'svg':
        charts.renderer = 'svg'
        make_snapshot(snapshot, charts.render(), save_path, delay=delay)
    else:
        make_snapshot(snapshot, charts.render(), save_path, delay=delay, pixel_ratio=pixel_ratio)