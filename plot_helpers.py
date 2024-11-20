import matplotlib.pyplot as plt
from matplotlib import font_manager


def initialize_mpl_style():

    plt.style.use('fivethirtyeight')
    plt.rcParams.update(
        {'axes.facecolor': 'w',
         'axes.edgecolor': 'w',
         'patch.edgecolor': 'w',
         'figure.facecolor': 'w',
         'savefig.edgecolor': 'w',
         'savefig.facecolor': 'w',
         'xtick.major.size': 0,
         'ytick.major.size': 0,
         'axes.axisbelow': True
         }
    )

    font_dir = "assets/font"
    font_files = font_manager.findSystemFonts(fontpaths=font_dir)
    for font_file in font_files:
        font_manager.fontManager.addfont(font_file)

    plt.rcParams["font.family"] = "Lato"