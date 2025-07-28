import gzip

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cartopy.io.shapereader as shpreader
from pyproj import Transformer
from pathlib import Path
import cmocean

from tawes import get_json_data, read_json_data
from plot_helpers import initialize_mpl_style

mm_levels = [0.5, 10, 20, 30, 40, 50, 60, 70, 80, 90,
             100, 110]
ticklabs = [str(x) for x in mm_levels]
lat_min, lat_max = (48.04, 48.40)
lon_min, lon_max = (16.10, 16.67)


def get_shapes():
    shpfile_countries = Path('assets/SHP_countries/LINIEN_050Line.shp')
    shpfile_provinces = Path('assets/SHP_provinces/LINIEN_050Line.shp')
    shpfile_districts = Path(
        'assets/SHP_districts/STATISTIK_AUSTRIA_POLBEZ_20230101.shp')
    shpfile_rivers = Path('assets/SHP_rivers/Routen.shp')

    uba_crs = ccrs.LambertConformal(
        central_longitude=13.3333, central_latitude=47.5, false_easting=400000,
        false_northing=400000, standard_parallels=(46, 49))

    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfile_countries).geometries(), ccrs.epsg(31259))
    provinces = cf.ShapelyFeature(
        shpreader.Reader(shpfile_provinces).geometries(), ccrs.epsg(31259))
    districts = cf.ShapelyFeature(
        shpreader.Reader(shpfile_districts).geometries(), ccrs.epsg(31287)
    )
    rivers = cf.ShapelyFeature(
        shpreader.Reader(shpfile_rivers).geometries(), uba_crs)

    return countries, provinces, districts, rivers



def plot_inca(ds, time0, time1, title, station_data=None, filename=None,
              ax=None, add_colorbar=True, x=(605000, 647500), y=(464000, 505000),
              cbar_kwargs=None):

    ds = ds.sel(x=slice(x[0], x[1]), y=slice(y[0], y[1]))
    ds_sum = ds.sel(time=slice(time0, time1)).sum(dim="time")
    c50, c75, c100 = ('limegreen', 'orange', 'red')
    #'#fc7e7e', '#ff4f4f', '#ff1919'

    initialize_mpl_style()

    norm = mpl.colors.BoundaryNorm(boundaries=mm_levels, ncolors=256)
    cmap = mpl.colormaps['cmo.rain']
    # Set the first color to white/transparent for very small values
    cmap.set_under('w', alpha=1)

    if ax is None:
        fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.epsg(31287)),
                               figsize=(9, 6))
    if add_colorbar:
        if cbar_kwargs is None:
            cbar_kwargs = {'extend': 'both', 'label': 'l/m²'}
    else:
        cbar_kwargs = None

    pm = ds_sum['RR'].plot.imshow(
        transform=ccrs.epsg(31287), ax=ax, norm=norm, interpolation='bilinear',
        cbar_kwargs=cbar_kwargs, add_colorbar=add_colorbar, levels=mm_levels,
        cmap=cmap, xlim=x, ylim=y)

    if add_colorbar is True:
        pm.colorbar.set_ticks(ticks=mm_levels, labels=ticklabs)
        pm.colorbar.outline.set_edgecolor('k')
        pm.colorbar.outline.set_linewidth(2)

    ds_sum['RR'].plot.contour(
        transform=ccrs.epsg(31287), ax=ax, levels=[50, 75, 100],
        linewidths=3, colors=[c50, c75, c100], zorder=11)

    for threshold, color in zip([50, 75, 100], [c50, c75, c100]):
        area = (ds_sum['RR'] > threshold).sum().values
        areatext = f"> {threshold} mm: {area} km²"
        if threshold == 50:
            text = ax.text(
                0, 0.01, areatext, color=color, fontweight='bold', zorder=15,
                bbox={'facecolor': 'white', 'pad': 1.5, 'zorder': 14},
                transform=ax.transAxes)
        else:
            text = ax.annotate(
                areatext, xycoords=text, xy=(0, 1.25), verticalalignment="bottom",
                color=color, fontweight='bold', zorder=15,
                bbox={'facecolor': 'w', 'pad': 1.5, 'zorder': 14})

    maxtext = (f"INCA {ds_sum['RR'].max().values:.0f} mm\n"
               f"TAWES {station_data['rr'].max():.0f} mm")
    text = ax.text(
        1, 0, maxtext,
        color='k', fontweight='bold', ha="right", va="bottom", zorder=15,
        bbox={'facecolor': 'w', 'pad': 0.5, 'zorder': 14},
        transform=ax.transAxes)

    ax.set_title(title)

    if station_data is not None:
        # plot TAWES locations
        plot_tawes(ax, station_data, pm.get_cmap(), pm.norm)

    # plot background shapes
    countries, provinces, districts, rivers = get_shapes()
    ax.add_feature(countries, linewidth=3, edgecolor='r',
                   facecolor='none', zorder=3)
    ax.add_feature(provinces, linewidth=2.5, edgecolor='k',
                   facecolor='none', zorder=2)
    ax.add_feature(rivers, linewidth=2, edgecolor='dodgerblue',
                   facecolor='none', zorder=1)
    ax.add_feature(districts, linewidth=1, edgecolor='k',
                   facecolor='none', zorder=1, alpha=0.5)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=200)
        return

    return pm


def plot_tawes(ax, station_data, cmap, norm):

    lon, lat, rr = (station_data['lon'], station_data['lat'], station_data['rr'])
    ax.scatter(lon, lat, s=70, c=rr, transform=ccrs.PlateCarree(), cmap=cmap,
               norm=norm, edgecolors='w', linewidth=1.6, zorder=12)
    return ax

# def plot_tawes_data(ax, station_data):
#
#     for station in station_data:
#         lon, lat = (station_data[station]['lon'], station_data[station]['lat'])
#         rr = np.round(station_data[station]['rr'], 1)
#         ax.scatter(lon, lat, color='k',transform=ccrs.PlateCarree())
#         ax.text(lon, lat+0.004, rr, ha='center', fontsize=12,
#                 transform=ccrs.PlateCarree())
#     return


def read_binary_format(inca_path, time0, time1):

    if time0 > dt.datetime(2012, 1, 1):
        coords = xr.open_dataset("data/georef.nc")
        data = coords.set_coords(["lat", "lon"])
    else:
        transformer = Transformer.from_crs("epsg:31287", "epsg:4326")

        x = np.arange(100*1000, 701*1000, 1000)
        y = np.arange(250*1000, 601*1000, 1000)

        xx, yy = np.meshgrid(x, y)
        lon, lat = transformer.transform(xx, yy)

        data = xr.Dataset(
            coords={
            "x": ("x", x),
            "y": ("y", y),
            "lon": (("y", "x"), lon),
            "lat": (("y", "x"), lat)})

    times = list()
    n = int((time1-time0).total_seconds()/60/15)

    rr = np.empty((n+1, data.y.size, data.x.size))
    i = 0

    time = time0
    while time <= time1:
        file = Path(inca_path, f"INCA_RR-{time:%H%M}.asc.gz")
        with gzip.open(file, 'rb') as f:
            x = f.read()
            x = np.fromstring(x, sep=' ')
            rr_temp = np.reshape(x, (data.y.size, data.x.size))
            rr[i, :, :] = rr_temp

            times.append(time)

        time += dt.timedelta(minutes=15)
        i += 1

    data.coords["time"] = pd.to_datetime(times)
    data["RR"] = (["time", "y", "x"], rr)
    data.attrs['freq'] = '15m'

    return data

def export_to_netcdf():
    start = dt.datetime(2023, 4, 1)
    end = dt.datetime(2023, 5, 1)
    current_day = start
    encode = {'RR': {'zlib': True, 'complevel': 9}}

    while current_day < end:
        inca_path = Path("data", f"{current_day:%Y}", f"{current_day:%Y%m%d}")
        print("processing", current_day)
        time0 = current_day
        time1 = current_day + dt.timedelta(hours=23, minutes=59)

        daily_data = read_binary_format(inca_path, time0, time1)
        daily_data.to_netcdf(Path("output", "INCA", f"{current_day:%Y%m%d}.nc"), encoding=encode)

        current_day += dt.timedelta(days=1)
    return


if __name__ == '__main__':

    print("nothing to do?")

    #export_to_netcdf()

    # event_starttimes = [
    #     dt.datetime(2003, 5, 13, 13, 0),
    #     dt.datetime(2004, 7, 1, 14, 0),
    #     dt.datetime(2007, 8, 9, 17, 0),
    #     dt.datetime(2008, 4, 22, 14, 0),
    #     dt.datetime(2008, 5, 18, 13, 0),
    #     dt.datetime(2010, 5, 13, 14, 0),
    #     dt.datetime(2011, 6, 8, 11, 0),
    #     dt.datetime(2011, 7, 28, 10, 0)
    #
    #     # dt.datetime(2014, 4, 29, 13, 0),
    #     # dt.datetime(2014, 5, 24, 13, 0),
    #     # dt.datetime(2018, 5, 2, 19, 0),
    #     # dt.datetime(2021, 7, 17, 18, 0),
    #     # dt.datetime(2024, 8, 17, 14, 0)
    # ]

    # for time0 in event_starttimes:
    #     time1 = time0 + dt.timedelta(hours=2)
    #
    #     if Path(f"data/{time0:%Y%m%d}/tawes.json").exists():
    #         pass
    #     else:
    #         get_json_data(time0, time1)
    #     station_data = read_json_data(f"data/{time0:%Y%m%d}/tawes.json")
    #     ds = read_binary_format(f"data/{time0:%Y%m%d}",
    #                             time0 + dt.timedelta(minutes=15),
    #                             time1)
    #     plot_title = f"{time0:%Y%m%d %H%M}-{time1:%H%M} UTC"
    #     plot_inca(ds, time0, time1, plot_title,
    #               filename=f'output/INCA_RR_{time0:%Y%m%d}.png',
    #               station_data=station_data)