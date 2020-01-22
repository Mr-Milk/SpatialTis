from skimage.external import tifffile

def tiff_header(img):
    with tifffile.TiffFile(img) as tif:
        for p in tif.pages:
            print(p.tags)

def diff_dicts_values(dict1, dict2):
    k1 = dict1.keys()

    for k in k1:
        diff1 = dict1[k] - dict2[k]
        diff2 = dict2[k] - dict1[k]
        if len(diff1) > 0:
            print(f'In dict1 not in dict2 {dff1}')
        if len(diff2) > 0:
            print(f'In dict2 not in dict1 {dff2}')

def plot_polygons(polygons):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in polygons:
        ax.plot(*i.exterior.xy, color='#6699cc', alpha=0.7,
            linewidth=3, solid_capstyle='round', zorder=2)
    ax.set_title('Polygons')