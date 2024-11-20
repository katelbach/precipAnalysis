import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature


import requests
import matplotlib.transforms as transforms

import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

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


def plot_radiosonde():

        if fig is None:
            fig = plt.figure(figsize=(9, 9))

        if time.minute != 0:
            raise ValueError('Radiosonde not available at minute: ',
                             time.minute)

        if time.hour not in [0, 6, 12, 18]:
            raise ValueError('Radiosonde not available at hour: ', time.hour)

        url = reader_radiosonde.construct_url(time)
        r = requests.get(url, allow_redirects=True)
        html = r.text
        df = reader_radiosonde.get_radiosonde_data_from_html(html)

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
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])

        # Calculate the parcel profile.
        parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')

        # Create a new figure. The dimensions here give a good aspect ratio
        skew = SkewT(fig, rotation=45, subplot=subplot)

        # Plot the data using normal plotting functions, in this case using
        # log scaling in Y, as dictated by the typical meteorological plot
        skew.plot(p, T, 'r')
        skew.plot(p, Td, 'g')
        skew.plot_barbs(p[p.m > 110][::80], u[p.m > 110][::80],
                        v[p.m > 110][::80])

        skew.ax.set_ylim(1000, 100)
        skew.ax.set_xlim(-30, 50)

        skew.ax.set_xlabel('Temperature [°C]')
        skew.ax.set_ylabel('Pressure [hPa]')

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

        # Create a hodograph
        # ax_hod = inset_axes(skew.ax, '30%', '30%', loc=1)
        ax_hod = skew.ax.inset_axes([0.65, 0.72, 0.26, 0.26])

        h = Hodograph(ax_hod, component_range=40.)
        h.add_grid(increment=20)
        h.plot_colormapped(u[p.m > 300][::30], v[p.m > 300][::30],
                           hght[p.m > 300][::30])

        for spine in ax_hod.spines.values():
            spine.set_linewidth(0.6)
            spine.set_color('black')

        for tick in [0, 20, 40]:
            ax_hod.text(tick - 4, -8, tick)

        ax_hod.set_xlim([-15, 45])

        ax_hod.get_xaxis().set_visible(False)
        ax_hod.get_yaxis().set_visible(False)

        if title:
            skew.ax.set_title('Vienna', loc='left', fontsize=14)
            skew.ax.set_title(f"{time:%d %b %Y %H%M UTC}", loc='right',
                              fontsize=14)
    return


if __name__ == '__main__':
    figure1()