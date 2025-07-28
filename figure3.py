# FIG 3: Estimation of a 100-yr return period event with CI
# 3 subplots in one row

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cmocean  # import for colormap
from pathlib import Path

from plot_helpers import initialize_mpl_style
from inca import get_shapes
from tawes import read_json_data


mm_levels = [0.5, 10, 20, 30, 40, 50, 60, 70, 80, 90,
             100, 110]
ticklabs = [str(x) for x in mm_levels]

def figure3():

    inca_ci = xr.open_dataset("data/review/inca_ci.nc")

    initialize_mpl_style()

    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=1, ncols=3, figsize=(14, 5),
        subplot_kw=dict(projection=ccrs.epsg(31287)))

    norm = mpl.colors.BoundaryNorm(boundaries=mm_levels, ncolors=256)
    cmap = mpl.colormaps['cmo.rain']
    # Set the first color to white/transparent for very small values
    cmap.set_under('w', alpha=1)

    cbar_kwargs = {'extend': 'both', 'label': '2 hour precipitation [mm]',
                   'drawedges': True}
    title_kwargs = {'color': 'k', 'fontweight': 'bold', 'fontsize': 17,
                    'ha': 'left', 'va': 'top', 'zorder': 15,
                    'bbox': {'facecolor': 'w', 'pad': 1, 'zorder': 14}}

    x = (605000, 645500)
    y = (465000, 500000)

    inca_ci.isel(CI=0)['prec'].plot.imshow(
        transform=ccrs.epsg(31287), ax=ax1, norm=norm, interpolation='bilinear',
        add_colorbar=False, levels=mm_levels,
        cmap=cmap, xlim=x, ylim=y)
    ax1.text(0, 1, "Lower CI", transform=ax1.transAxes, **title_kwargs)

    inca_ci.isel(CI=1)['prec'].plot.imshow(
        transform=ccrs.epsg(31287), ax=ax2, norm=norm, interpolation='bilinear',
        add_colorbar=False, levels=mm_levels, cmap=cmap, xlim=x, ylim=y)
    ax2.text(0, 1, "Mean", transform=ax2.transAxes, **title_kwargs)

    im3 = inca_ci.isel(CI=2)['prec'].plot.imshow(
        transform=ccrs.epsg(31287), ax=ax3, norm=norm, interpolation='bilinear',
        add_colorbar=False, levels=mm_levels, cmap=cmap, xlim=x, ylim=y)
    ax3.text(0, 1, "Upper CI", transform=ax3.transAxes, **title_kwargs)

    plt.subplots_adjust(left=0.02, bottom=0.05, top=0.95, wspace=0.1, right=0.9)
    cbar_ax = fig.add_axes((0.91, 0.15, 0.02, 0.7))
    cbar = fig.colorbar(im3, cax=cbar_ax, **cbar_kwargs)

    cbar.set_ticks(ticks=mm_levels, labels=ticklabs)
    cbar.outline.set_edgecolor('k')
    cbar.outline.set_linewidth(2)

    # add shapes
    countries, provinces, districts, rivers = get_shapes()
    for ax in [ax1, ax2, ax3]:
        ax.add_feature(countries, linewidth=3, edgecolor='r',
                       facecolor='none', zorder=3)
        ax.add_feature(provinces, linewidth=2.5, edgecolor='k',
                       facecolor='none', zorder=2)
        ax.add_feature(rivers, linewidth=2, edgecolor='dodgerblue',
                       facecolor='none', zorder=1)
        ax.add_feature(districts, linewidth=1, edgecolor='k',
                       facecolor='none', zorder=1, alpha=0.5)

    # unset auto-generated title
    for ax in [ax1, ax2, ax3]:
        ax.set_title("")

    # plot TAWES locations
    station_data = read_json_data(f"data/20240817/tawes.json")
    station_data = pd.DataFrame.from_dict(station_data, orient='index')

    lon, lat = (station_data['lon'], station_data['lat'])
    for ax in [ax1, ax2, ax3]:
        ax.scatter(lon, lat, s=40, c='w', transform=ccrs.PlateCarree(),
                   edgecolors='k', linewidth=1.6, zorder=12)

    Path("output").mkdir(exist_ok=True)
    plt.savefig("output/figure3.png", dpi=300, bbox_inches='tight')
    plt.close()

    return

if __name__ == "__main__":
    figure3()