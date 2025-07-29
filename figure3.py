# FIG 3: Estimation of a 100-yr return period event with CI
# 3 subplots in one row

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cmocean  # import for colormap
from pathlib import Path

from metpy.io import station_data

from plot_helpers import initialize_mpl_style
from inca import get_shapes


mm_levels = [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
ticklabs = [str(x) for x in mm_levels]

station_metadata = {
    'Wien Hohe Warte': (48.24861,16.35639),
    'Wien Mariabrunn': (48.20694351, 16.2294445),
    'Wien Unterlaa': (48.125, 16.41944313),
    'Wien Jubiläumswarte': (48.2211113, 16.26527786),
    'Brunn am Gebirge': (48.10694504, 16.26999855),
    'Wien Innere Stadt': (48.19833374, 16.36694336),
    'Wien Donaufeld': (48.25722122, 16.43138885),
    'Schwechat Flughafen': (48.11750031, 16.58138847),
    'Groß Enzersdorf': (48.19972229, 16.55916595),
    'Wien Stammersdorf': (48.30583191, 16.40555573),
    'Langenlebarn': (48.32388687, 16.11805534)
}

def figure3():

    inca_ci = xr.open_dataset("data/review/inca_ci.nc")
    station_ci = pd.read_csv("data/review/stations_ci_100year.csv")

    station_ci[['lat', 'lon']] = station_ci['station'].map(
        station_metadata).apply(pd.Series)

    initialize_mpl_style()

    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=1, ncols=3, figsize=(14, 5),
        subplot_kw=dict(projection=ccrs.epsg(31287)))

    norm = mpl.colors.BoundaryNorm(boundaries=mm_levels, ncolors=256,
                                   extend='both')
    cmap = mpl.colormaps['viridis_r']

    cbar_kwargs = {'extend': 'both', 'label': '2-hour precipitation [mm]',
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
    ax1.text(0, 1, "Lower CI (2.5 %)", transform=ax1.transAxes, **title_kwargs)

    inca_ci.isel(CI=1)['prec'].plot.imshow(
        transform=ccrs.epsg(31287), ax=ax2, norm=norm, interpolation='bilinear',
        add_colorbar=False, levels=mm_levels, cmap=cmap, xlim=x, ylim=y)
    ax2.text(0, 1, "Mean", transform=ax2.transAxes, **title_kwargs)

    im3 = inca_ci.isel(CI=2)['prec'].plot.imshow(
        transform=ccrs.epsg(31287), ax=ax3, norm=norm, interpolation='bilinear',
        add_colorbar=False, levels=mm_levels, cmap=cmap, xlim=x, ylim=y)
    ax3.text(0, 1, "Upper CI (97.5 %)", transform=ax3.transAxes, **title_kwargs)

    for i, ax in enumerate([ax1, ax2, ax3]):
        axtext = (f"Max {inca_ci.isel(CI=i)['prec'].max().values:.0f} mm\n"
                   f"Avg {inca_ci.isel(CI=i)['prec'].mean():.0f} mm\n"
                   f"Min  {inca_ci.isel(CI=i)['prec'].min():.0f} mm")
        ax.text(
            0.0, 0, axtext,
            color='k', fontweight='bold', ha="left", va="bottom", zorder=15,
            fontsize=13,
            bbox={'facecolor': 'w', 'pad': 0.5, 'zorder': 14, 'alpha': 1.0},
            transform=ax.transAxes)

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
        ax.add_feature(districts, linewidth=1, edgecolor='k',
                       facecolor='none', zorder=1, alpha=0.5)
        ax.add_feature(rivers, linewidth=2, edgecolor='dodgerblue',
                       facecolor='none', zorder=3)

    # unset auto-generated title
    for ax in [ax1, ax2, ax3]:
        ax.set_title("")

    # plot TAWES locations
    lon, lat, mean, c_low, c_high = (
        station_ci['lon'], station_ci['lat'],
        station_ci['mean_estimate'], station_ci['bound.0.025'],
        station_ci['bound.0.975'])

    for ax, val in zip([ax1, ax2, ax3], [c_low, mean, c_high]):
        ax.scatter(lon, lat, s=70, c=val, transform=ccrs.PlateCarree(),
                   cmap=cmap, norm=norm, edgecolors='w', linewidth=1.6,
                   zorder=12)

    Path("output").mkdir(exist_ok=True)
    plt.savefig("output/figure3.png", dpi=300, bbox_inches='tight')
    plt.close()

    return

if __name__ == "__main__":
    figure3()