import datetime as dt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

from tawes import read_csv_data, evt_analysis
from plot_helpers import initialize_mpl_style

# FIG1: 2-hour data Wien H.W.
# a) annual maxima of 2h precipitation
# b) return value estimate -> annual maxima with GEV as distribution

def figure1():

    obs_timeseries = read_csv_data()

    obs_quantiles = pd.read_csv(
        "data/review/vienna_hohe-warte_emp_quantiles.csv")
    model_quantiles = pd.read_csv(
        "data/review/vienna_hohe-warte_theo_quantiles_and_ci.csv")

    # FIG 1: Return period analysis for Vienna H.W.
    initialize_mpl_style()
    fontsize_labels = 12
    xlim = [1, 800]
    ylim = [0, 120]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

    # SUBPLOT A: Time series of daily 2hr RR maximum
    rr_2h = obs_timeseries['2H_SUM'].dropna()
    ax1.plot(rr_2h.index, rr_2h, ls="-", color="cornflowerblue",
             lw=0.25, zorder=10)

    min_year = rr_2h.index.year.min()
    max_year = rr_2h.index.year.max()

    # Configure axes
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(color='k', lw=0.3, alpha=0.5, zorder=5)

    # plot vertical lines every 10 years
    for yr in range(int(np.ceil(min_year/10) * 10), max_year+1, 10):
        ax1.axvline(dt.datetime(yr, 1, 1), lw=0.3, color="k", alpha=0.5,
                   zorder=5)

    ax1.set_xlim(dt.datetime(min_year, 1, 1), dt.datetime(max_year+1, 12, 1))
    ax1.set_ylim(ylim)

    ax1.set_xlabel('Year', fontsize=fontsize_labels)
    ax1.set_ylabel('mm', fontsize=fontsize_labels)
    ax1.text(0, 1, '(a) Max 2h RR', ha='left', va='top',
             transform=ax1.transAxes,
             zorder=12, fontweight='bold',
             bbox=dict(facecolor='w', edgecolor='none'))


    # SUBPLOT B: Return periods for Hohe Warte
    ax2.set_xscale('log')
    ax2.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:,.0f}"))

    ax2.fill_between(model_quantiles["rp"], model_quantiles['bound.0.25'],
                     model_quantiles['bound.0.75'], alpha=.4,
                     color='cornflowerblue')
    ax2.fill_between(model_quantiles["rp"], model_quantiles['bound.0.025'],
                     model_quantiles['bound.0.975'], alpha=.2,
                     color='cornflowerblue')
    ax2.plot(model_quantiles["rp"], model_quantiles['bound.0.025'], "--",
             color='cornflowerblue', linewidth=1)
    ax2.plot(model_quantiles["rp"], model_quantiles['bound.0.975'], "--",
             color='cornflowerblue', linewidth=1)
    ax2.plot(model_quantiles['rp'], model_quantiles['mean_estimate'],
             color='r', linewidth=2)

    ax2.axhline(107, linestyle='--', color='k', linewidth=2.5)
    ax2.scatter(obs_quantiles['rp'], obs_quantiles['obs'], s=15, color='k',
                edgecolor='w')

    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim)

    ax2.set_xticks([1, 2, 5, 10, 30, 50, 100, 200, 500])
    ax2.set_xlabel('Years', fontsize=fontsize_labels)
    ax2.set_ylabel(None)
    ax2.text(0, 1, '(b) Return periods 2h RR', ha='left', va='top',
             transform=ax2.transAxes,
             zorder=12, fontweight='bold',
             bbox=dict(facecolor='w', edgecolor='none'))

    plt.tight_layout()

    Path("output").mkdir(exist_ok=True)
    plt.savefig("output/figure1.pdf", format="pdf", dpi=300,
                bbox_inches='tight')
    plt.close()

    return


def plot_timeseries(df, ax=None):

    rr_2h = (
        df['2H_SUM']
        .sort_index(ascending=True)
        .astype(float)
        .dropna())

    rr_2h_annualmax = rr_2h.groupby(rr_2h.index.year).max()
    idxmax = rr_2h.groupby(rr_2h.index.year).idxmax()
    rr_2h_annualmax.index = idxmax

    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

    # Configure axes
    ax.xaxis.grid(False)
    ax.yaxis.grid(color='k', lw=0.3, alpha=0.5, zorder=5)

    # Plot signal time series
    ax.plot(rr_2h.index, rr_2h, ls="-", color="cornflowerblue",
            lw=0.25, zorder=10)

    # Plot extreme events
    ax.scatter(rr_2h_annualmax.index, rr_2h_annualmax.values, s=25, lw=0.5,
               edgecolor="w", facecolor="r", zorder=15)

    min_year = rr_2h.index.year.min()
    max_year = rr_2h.index.year.max()

    # plot vertical lines every 10 years
    for yr in range(int(np.ceil(min_year/10) * 10), max_year+1, 10):
        ax.axvline(dt.datetime(yr, 1, 1), lw=0.3, color="k", alpha=0.5,
                   zorder=5)

    ax.set_xlim(dt.datetime(min_year, 1, 1), dt.datetime(max_year+1, 12, 1))

    return ax


def plot_return_levels(model, ax=None):

    return_period_steps = [1.01, 1.05, 1.1, 1.2, 1.4, 1.55, 1.68, 2, 3, 4, 5,
                           7, 10, 15, 25, 35, 50, 75, 100, 150, 200, 300, 400,
                           500, 600, 700]

    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)

    model.plot_return_values(ax=ax, alpha=0.95,
                             return_period=return_period_steps)
    ax.set_xticks([1, 2, 5, 10, 30, 50, 100, 200, 500])
    ax.xaxis.grid(False, which='both')
    ax.xaxis.grid()

    return ax

if __name__ == '__main__':
    figure1()
