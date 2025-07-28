# FIG 3: Estimation of a 100-yr return period event with CI
# 3 subplots in one row

import datetime as dt
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs


mm_levels = [0.1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90,
             100, 110]

def figure3():

    inca_ci = xr.open_dataset("data/review/inca_ci.nc")

    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=1, ncols=3, figsize=(12, 5),
        subplot_kw=dict(projection=ccrs.epsg(31287)))

    norm = colors.BoundaryNorm(boundaries=mm_levels, ncolors=256)

    inca_ci.isel(ci=0).plot.imshow(
        transform=ccrs.epsg(31287), ax=ax1, norm=norm, interpolation='bilinear',
        cbar_kwargs=cbar_kwargs, add_colorbar=True, levels=mm_levels,
        cmap='cmo.rain', xlim=x, ylim=y)



    return

if __name__ == "__main__":
    figure3()