import xarray as xr
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units

from plot_helpers import initialize_mpl_style
from radiosonde import get_profile


def figure1():
    filename_geopot = "data/20240817/geopot.grib"
    geopot = xr.open_dataset(filename_geopot, engine='cfgrib')
    filename_mslp = "data/20240817/mslp.grib"
    mslp = xr.open_dataset(filename_mslp, engine='cfgrib')

    geopot = prepare_grib_dataset(geopot, (-50, 65), (70, 20))
    mslp = prepare_grib_dataset(mslp, (-50, 65), (70, 20))

    geopot500 = geopot['z'].sel(isobaricInhPa=500) / 100
    mslp = mslp['msl'] / 100
    aeqrtp = (geopot['z'].sel(isobaricInhPa=500) -
              geopot['z'].sel(isobaricInhPa=1000)) / 100

    initialize_mpl_style()

    fig = plt.figure(figsize=(14, 6))
    # hack for additional wspace between colorbar ax and radiosonde plot ax
    gs = fig.add_gridspec(1, 4, width_ratios=[1.2, 0.05, 0.13, 1])
    projection = ccrs.LambertConformal(
        central_longitude=15.0, central_latitude=46.0)

    # Synoptic overview
    ax1 = fig.add_subplot(gs[0], projection=projection)

    pm = aeqrtp.plot.contourf(
        ax=ax1, cmap='gist_rainbow_r', levels=np.arange(480, 600, 4),
        transform=ccrs.PlateCarree(), add_colorbar=False)

    # add colorbar
    ax_cbar = fig.add_subplot(gs[1])
    cbar = fig.colorbar(pm, cax=ax_cbar, extend='both', shrink=0.5)
    cbar.set_label('AeqRTP H500-H1000 [gpdm]', labelpad=-65)

    # add 500 hPa geopot contours
    con1 = geopot500.plot.contour(
        ax=ax1, levels=np.arange(400, 700, 8), linewidths=2,
        colors='k', transform=ccrs.PlateCarree())
    ax1.clabel(con1, fontsize=14, inline=1, inline_spacing=-8, fmt='%i',
              rightside_up=True)

    # add marker at the location of Vienna
    ax1.scatter(16.36, 48.25, marker='o', s=60, color='red',
             transform=ccrs.PlateCarree())

    # add MSLP contours
    con2 = ax1.contour(
        mslp.longitude, mslp.latitude, gaussian_filter(mslp, sigma=2),
        levels = np.arange(950, 1050, 5), linewidths=2,
        colors='w', transform=ccrs.PlateCarree())
    ax1.clabel(con2, con2.levels, fontsize=16, inline=1, inline_spacing=-13,
               fmt='%i', rightside_up=True)

    ax1.text(0, 1, f"a) RelTop/MSLP/500 hPa Geopot", fontweight='bold',
             transform=ax1.transAxes, verticalalignment='top', fontsize=20,
             horizontalalignment='left', bbox={'facecolor': 'w', 'pad': 1},
             zorder=11)
    ax1.set_title("")

    ax1.set_extent([-20, 40, 67, 27], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.75)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.5)

    radiosonde_time = dt.datetime(2024, 8, 17, 12, 0)
    ax2 = plot_radiosonde(radiosonde_time, fig=fig, subplot=gs[3])

    ax2.set_ylim(1000, 105)
    ax2.set_xlim(-30, 50)

    ax2.set_xlabel('Temperature [°C]')
    ax2.set_ylabel('Pressure [hPa]', labelpad=-5)

    ax2.text(0, 1, f"b) Vienna radiosonde", fontweight='bold',
             transform=ax2.transAxes, verticalalignment='top', fontsize=20,
             horizontalalignment='left', bbox={'facecolor': 'w', 'pad': 1},
             zorder=11)


    plt.subplots_adjust(wspace=0.1, top=0.99, bottom=0.1, left=0.01,
                        right=0.97)
    # debugging only
    # plt.show()

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


def plot_radiosonde(time, fig=None, subplot=None):

    if fig is None:
        fig = plt.figure(figsize=(8, 6))
    if subplot is None:
        subplot = 111

    df = get_profile(time)

    # Drop any rows with all NaN values for T, Td, winds
    df = df.dropna(subset=('TEMP', 'DWPT', 'DRCT', 'SPED'),
                   how='all').reset_index(drop=True)

    # filter for same pressure values
    p = df['PRES'].values * units.hPa
    p_diff = np.diff(p, append=p[-1])
    mask = p_diff != 0
    df = df[mask]

    p = df['PRES'].values * units.hPa
    T = df['TEMP'].values * units.degC
    Td = df['DWPT'].values * units.degC
    hght = df['HGHT'].values * units.m
    wind_speed = df['SPED'].values * units('m/s')
    wind_speed = wind_speed.to('knots')
    wind_dir = df['DRCT'].values * units.degrees
    u, v = mpcalc.wind_components(wind_speed, wind_dir)

    # Calculate the LCL
    # lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])

    # Calculate the parcel profile.
    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')

    # Create a new figure. The dimensions here give a good aspect ratio
    skew = SkewT(fig, rotation=45, subplot=subplot)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.plot(p, T, 'r', linewidth=2.5)
    skew.plot(p, Td, 'g', linewidth=2.5)
    skew.plot_barbs(p[p.m > 110][::80], u[p.m > 110][::80],
                    v[p.m > 110][::80])

    # Plot LCL temperature as black dot
    # skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

    # Plot the parcel profile as a black line
    skew.plot(p, parcel_prof, 'k', linewidth=1.5, alpha=0.8)

    # Shade areas of CAPE and CIN
    # skew.shade_cin(p, T, parcel_prof, Td)
    skew.shade_cape(p, T, parcel_prof, alpha=0.2)

    # Plot a zero degree isotherm
    # skew.ax.axvline(0, color='c', linestyle='--', linewidth=1.5)

    # Add the relevant special lines
    skew.plot_dry_adiabats(alpha=0.4, linewidth=1)
    skew.plot_moist_adiabats(alpha=0.4, linewidth=1)
    skew.plot_mixing_lines(alpha=0.4, linewidth=1)

    # add height in [m] for certain pressure levels
    critical_pressure_levels = np.array([900, 850, 700, 500, 400, 300, 200])
    pressure_idx = _closest_argmin(critical_pressure_levels, p.m)
    p_critical = p.m[pressure_idx]
    hght_critical = hght.m[pressure_idx]

    trans = transforms.blended_transform_factory(
        skew.ax.transAxes, skew.ax.transData)

    for pcrit, hcrit in zip(p_critical, hght_critical):
        skew.ax.text(0.02, pcrit, str(hcrit), color='grey', transform=trans,
                     fontsize=8)

    return skew.ax


def _closest_argmin(A, B):
    L = B.size
    sidx_B = B.argsort()
    sorted_B = B[sidx_B]
    sorted_idx = np.searchsorted(sorted_B, A)
    sorted_idx[sorted_idx==L] = L-1
    mask = (sorted_idx > 0) & \
    ((np.abs(A - sorted_B[sorted_idx-1]) < np.abs(A - sorted_B[sorted_idx])) )
    return sidx_B[sorted_idx-mask]


if __name__ == '__main__':
    figure1()