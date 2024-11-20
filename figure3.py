from pathlib import Path
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from inca import plot_inca, get_json_data, read_binary_format
from tawes import read_json_data

# FIG2: 5 subplots, left big: Main event, right small: four next-largest events (sorted by max 2h RR)

def figure3():

    # main event
    time0 = dt.datetime(2024, 8, 17, 14, 0)
    time1 = time0 + dt.timedelta(hours=2)
    inca_ds = get_data(time0, time1)

    station_data = read_json_data(f"data/{time0:%Y%m%d}/tawes.json")
    station_data = pd.DataFrame.from_dict(station_data, orient='index')

    fig = plt.figure(layout='constrained', figsize=(12, 6))
    subfigs = fig.subfigures(1, 2, wspace=0.07, width_ratios=[1, 1.2])

    ax1 = subfigs[0].subplots(1, 1, subplot_kw=dict(projection=ccrs.epsg(31287)))
    pm = plot_inca(
        inca_ds, time0, time1, "", filename=None, ax=ax1, add_colorbar=True,
        cbar_kwargs={'extend': 'both', 'label': '2 hour precipitation [mm]',
                     'location': 'bottom'})
    plot_tawes_locations(ax1, station_data, pm.get_cmap(), pm.norm)
    ax1.text(0.01, 1, f"a) {time0:%Y-%m-%d %H:%M}", fontweight='bold',
            transform=ax1.transAxes, verticalalignment='top', fontsize=18,
            horizontalalignment='left', bbox={'facecolor': 'w', 'pad': 1})
    ax1.axis('off')

    old_events = [
        dt.datetime(2014, 5, 24, 13, 0), # 62,1 mm
        dt.datetime(2021, 7, 17, 18, 0), # 57,8 mm
        dt.datetime(2010, 5, 13, 14, 0), # 55,2 mm
        dt.datetime(2008, 5, 18, 13, 0) # 51,6 mm
    ]
    # next-ranking smaller events (same colorbar!)
    axs2 = subfigs[1].subplots(2, 2,
                               subplot_kw=dict(projection=ccrs.epsg(31287)))

    label = ["b)", "c)", "d)", "e)"]
    for i, ax in enumerate(axs2.flatten()):
        time0 = old_events[i]
        time1 = time0 + dt.timedelta(hours=2)
        inca_ds = get_data(time0, time1)
        station_data = read_json_data(f"data/{time0:%Y%m%d}/tawes.json")
        station_data = pd.DataFrame.from_dict(station_data, orient='index')

        pm = plot_inca(inca_ds, time0, time1, "", filename=None, ax=ax,
                  add_colorbar=False, y=(460000, 503000))
        plot_tawes_locations(ax, station_data, pm.get_cmap(), pm.norm)
        ax.text(0.01, 1, f"{label[i]} {time0:%Y-%m-%d %H:%M}", fontweight='bold',
                transform=ax.transAxes, verticalalignment='top', fontsize=18,
                horizontalalignment='left', bbox={'facecolor': 'w', 'pad': 1})

    plt.savefig("output/figure3.png", dpi=200, bbox_inches='tight')
    plt.close()

    return

def plot_tawes_locations(ax, station_data, cmap, norm):

    lon, lat, rr = (station_data['lon'], station_data['lat'], station_data['rr'])
    ax.scatter(lon, lat, s=70, c=rr, transform=ccrs.PlateCarree(), cmap=cmap,
               norm=norm, edgecolors='w', linewidth=1.6)
    return ax


def get_data(time0, time1):
    ds = read_binary_format(
        f"data/{time0:%Y%m%d}", time0 + dt.timedelta(minutes=15), time1)
    return ds


if __name__ == '__main__':
    figure3()
