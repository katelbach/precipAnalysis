import matplotlib.pyplot as plt


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
    #plt.rcParams["font.family"] = "Arial"