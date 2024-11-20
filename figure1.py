import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from plot_helpers import initialize_mpl_style


def figure1():
    filename_geopot = "data/20240817/geopot.grib"
    geopot = xr.open_dataset(filename_geopot, engine='cfgrib')
    filename_mslp = "data/20240817/mslp.grib"
    mslp = xr.open_dataset(filename_mslp, engine='cfgrib')

    geopot = prepare_grib_dataset(geopot, (-50, 55), (70, 20))
    mslp = prepare_grib_dataset(mslp, (-50, 55), (70, 20))

    geopot500 = geopot['z'].sel(isobaricInhPa=500) / 100
    mslp = mslp['msl'] / 100
    aeqrtp = (geopot['z'].sel(isobaricInhPa=500) -
              geopot['z'].sel(isobaricInhPa=1000)) / 100

    initialize_mpl_style()

    fig = plt.figure(figsize=(13, 6), constrained_layout=True)
    projection = ccrs.LambertConformal(
        central_longitude=15.0, central_latitude=46.0)

    # Synoptic overview
    ax1 = fig.add_subplot(121, projection=projection)

    pm = aeqrtp.plot.contourf(
        ax=ax1, cmap='gist_rainbow_r', levels=np.arange(480, 600, 4),
        transform=ccrs.PlateCarree(), add_colorbar=True,
        cbar_kwargs={'location': 'bottom', 'extend': 'both', 'shrink': 0.7,
                     'label': 'AeqRTP H500-H1000 [gpdm]'})
    con1 = geopot500.plot.contour(
        ax=ax1, levels=np.arange(400, 700, 8), linewidths=2,
        colors='k', transform=ccrs.PlateCarree())
    ax1.clabel(con1, fontsize=14, inline=1, inline_spacing=-8, fmt='%i',
              rightside_up=True)

    con2 = ax1.contour(
        mslp.longitude, mslp.latitude, gaussian_filter(mslp, sigma=2),
        levels = np.arange(950, 1050, 5), linewidths=2,
        colors='w', transform=ccrs.PlateCarree())
    ax1.clabel(con2, con2.levels, fontsize=16, inline=1, inline_spacing=-13,
               fmt='%i', rightside_up=True)

    ax1.text(0, 1, f"a) RelTop/MSLP/500 hPa Geopot", fontweight='bold',
             transform=ax1.transAxes, verticalalignment='top', fontsize=18,
             horizontalalignment='left', bbox={'facecolor': 'w', 'pad': 1},
             zorder=11)
    ax1.set_title("")

    ax1.set_extent([-20, 40, 60, 30], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.75)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.5)

    # CAPE + PWAT
    filename_conv = "data/20240817/pwat_cape.grib"
    conv = xr.open_dataset(filename_conv, engine='cfgrib')
    conv = prepare_grib_dataset(conv, (-50, 55), (70, 20))

    ax2 = fig.add_subplot(122, projection=projection)
    cape_levels = [0, 100, 200, 300, 500, 750, 1000, 1250, 1500, 2000, 2500]
    conv['cape'].plot.contourf(
        ax=ax2, levels=cape_levels, cmap='jet', add_colorbar=True,
        cbar_kwargs={'location': 'bottom', 'extend': 'both', 'shrink': 0.7,
                     'label': 'CAPE [J/kg]'},
        transform=ccrs.PlateCarree())
    con = ax2.contourf(
        conv.longitude, conv.latitude, gaussian_filter(conv['tcwv'], sigma=1.5),
        levels=[40, 90], hatches=['//'], colors='none',
        transform=ccrs.PlateCarree())
    con.set_edgecolor('w')

    con2 = ax2.contour(
        conv.longitude, conv.latitude, gaussian_filter(conv['tcwv'], sigma=1.5),
        levels=[40], colors='w', transform=ccrs.PlateCarree(),
        linewidths=2)
    plt.rcParams.update({'hatch.color': 'w'})

    ax2.set_extent([0, 30, 55, 40], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.75)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax2.set_title("")
    ax2.text(0, 1, f"b) CAPE/PWAT", fontweight='bold',
             transform=ax2.transAxes, verticalalignment='top', fontsize=18,
             horizontalalignment='left', bbox={'facecolor': 'w', 'pad': 1},
             zorder=11)

    # debugging only
    #plt.show()

    plt.savefig("output/figure1.png", dpi=200, bbox_inches='tight')
    plt.close()

    return

def prepare_grib_dataset(ds, xlim, ylim):
    ds['longitude'] = xr.where(
        ds.longitude > 180, ds.longitude-360, ds.longitude)

    ds = ds.roll(longitude=int(len(ds['longitude']) / 2) - 1,
                 roll_coords=True)
    ds = ds.sel(longitude=slice(xlim[0], xlim[1]),
                latitude=slice(ylim[0], ylim[1]))

    return ds


if __name__ == '__main__':
    figure1()