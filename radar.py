import datetime as dt
from pathlib import Path
import matplotlib.pyplot as plt
import xradar as xd
import cartopy.crs as ccrs
import pyart

from pyraps.io import odim_h5
from inca import get_shapes


time0 = dt.datetime(2014, 5, 24, 13, 0)
time1 = dt.datetime(2014, 5, 24, 15, 0)

radardir = Path("data/20140524")

# vcp_list = list()
time = time0
countries, provinces, districts, rivers = get_shapes()

while time <= time1:

    filename = f"PARA01_ACGR_{time:%Y%m%d%H%M}.hdf"
    vcp = odim_h5.read_odim_h5(Path(radardir, filename))

    vcp_ds = vcp.pyraps.datatree2set(dict(
        tolerance=0.25, method='nearest', start_angle=0, stop_angle=360,
        angle_res=1, direction=1))
    vcp_ds = vcp_ds.xradar.georeference()
    vcp_ds.coords["starttime"] = time

    vcp_ds = vcp_ds.where((vcp_ds['DBZH'] > 0) & (vcp_ds['DBZH'] < 80))
    # vcp_list.append(vcp_ds)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection=ccrs.epsg(31287))
    proj_crs = xd.georeference.get_crs(vcp_ds)
    cart_crs = ccrs.Projection(proj_crs)

    vcp_ds.isel(fixed_angle=-2)['DBZH'].plot(x="x", y="y", ax=ax, vmin=0,
                                             vmax=70, cmap='pyart_NWSRef',
                                             transform=cart_crs)
    ax.gridlines(draw_labels=True)
    ax.add_feature(countries, linewidth=2, edgecolor='r',
                   facecolor='none', zorder=3)
    ax.add_feature(provinces, linewidth=1, edgecolor='k',
                   facecolor='none', zorder=2)
    ax.add_feature(rivers, linewidth=2, edgecolor='dodgerblue',
                   facecolor='none', zorder=1)
    ax.add_feature(districts, linewidth=1, edgecolor='k',
                   facecolor='none', zorder=1, alpha=0.5)
    ax.set_extent([15.5, 16.8, 47.8, 48.7], crs=ccrs.PlateCarree())
    plt.savefig(f"output/RAU_20140524/DBZH_SWP2_{time:%Y%m%d%H%M}.png",
                bbox_inches='tight')
    plt.close()

    time += dt.timedelta(minutes=5)

# vcps = xr.concat(vcp_list, dim="starttime")
