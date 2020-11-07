paper_context = {
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.linewidth': 1.0,
    'grid.linewidth': 0.8,
    'lines.linewidth': 0.8,
    'lines.markersize': 3,
    'patch.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.minor.size': 2,
    'ytick.minor.size': 2
}


paper_context2 = {
    # Font
    # ==================================================================
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "axes.titleweight": "normal",
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.title_fontsize": 7,
    # Spines, grid, ticks
    # ==================================================================
    "axes.spines.bottom": True,
    "axes.spines.left": True,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 1.0,
    "grid.linewidth": 1.0,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 2,
    "ytick.minor.size": 2,
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "xtick.minor.pad": 2,
    "ytick.minor.pad": 2,
    "axes.labelpad": 4.0,
    "axes.titlepad": 4.0,
    "axes.xmargin": 0,
    "axes.ymargin": 0,
    # plotted lines, markers, patches
    # ==================================================================
    "lines.linewidth": 1,
    "patch.linewidth": 1,
    # it appears that the default marker style and size in ggplot and plotnine is hardcoded
    # https://plotnine.readthedocs.io/en/stable/generated/plotnine.geoms.geom_point.html#plotnine.geoms.geom_point
    # https://stackoverflow.com/questions/34638902/point-size-in-ggplot-2-0-0
    "lines.marker": "None",
    "lines.markeredgewidth": 0.5,
    # these should not be set for mpl, because they will currently interfere with
    # labeling the tick marks, eg for colorbars, see also:
    # https://github.com/matplotlib/matplotlib/issues/14546
    # 'lines.markeredgecolor': 'black',
    # 'lines.markerfacecolor': 'black',
    "lines.markersize": ((1.5 + 0.5) ** 2),
    "scatter.marker": "o",
    # legend
    # ==================================================================
    # in font size units
    # details: https://matplotlib.org/3.1.1/api/legend_api.html
    "legend.handlelength": 0.5,
    "legend.handleheight": 0.5,
    # vertical space between legend entries
    "legend.labelspacing": 0.5,
    # pad between axes and legend border
    "legend.borderaxespad": 0,
    # fractional whitspace inside legend border, still measured in font size units
    "legend.borderpad": 0,
    "legend.frameon": False,
    # constrained_layout
    # ==================================================================
    "figure.constrained_layout.use": False,
    "figure.constrained_layout.h_pad": 0,
    "figure.constrained_layout.w_pad": 0,
    "figure.constrained_layout.hspace": 0.05,
    "figure.constrained_layout.wspace": 0.05,
    # Export options
    # ==================================================================
    "figure.facecolor": "white",
    # In case you do need raster graphics, ProPlot sets the default rc[‘savefig.dpi’] to 1200 dots per inch, which is recommended by most journals as the minimum resolution for rasterized figures containing lines and text. See the configuration section for how to change any of these settings.
    # https://proplot.readthedocs.io/en/latest/basics.html
    "savefig.dpi": 1200,
    "figure.dpi": 180,
    # avoid exporting text as one path per letter in svg and pdf
    # - https://stackoverflow.com/questions/34387893/output-matplotlib-figure-to-svg-with-text-as-text-not-curves
    # - https://matplotlib.org/tutorials/introductory/customizing.html
    #   svg.fonttype: path      # How to handle SVG fonts:
    #       path: Embed characters as paths -- supported
    #            by most SVG renderers
    #       None: Assume fonts are installed on the
    #            machine where the SVG will be viewed.
    "svg.fonttype": "none",
    # Add "pdf.use14corefonts: True" in your configuration file to use only
    # the 14 PDF core fonts. These fonts do not need to be embedded; every
    # PDF viewing application is required to have them. This results in very
    # light PDF files you can use directly in LaTeX or ConTeXt documents
    # generated with pdfTeX, without any conversion.
    # These fonts are: Helvetica, Helvetica-Bold, Helvetica-Oblique,
    # Helvetica-BoldOblique, Courier, Courier-Bold, Courier-Oblique,
    # Courier-BoldOblique, Times-Roman, Times-Bold, Times-Italic,
    # Times-BoldItalic, Symbol, ZapfDingbats.
    # "pdf.use14corefonts": True,
}
