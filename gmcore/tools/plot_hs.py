#!/usr/bin/env python3

import matplotlib.pyplot as plt
import xarray as xr

f = xr.open_dataset('~/hs.128x72x26.h0.nc', decode_times=False)

ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

u = f.u.isel(time=slice(-1000, None)).mean(dim='time')
u = u.mean(dim='lon')

u.plot.contourf(ax=ax1, y='lev', yincrease=False, cmap='jet', levels=20)

t = f.t.isel(time=slice(-1000, None)).mean(dim='time')
t = t.mean(dim='lon')

t.plot.contourf(ax=ax2, y='lev', yincrease=False, cmap='jet', levels=20)

plt.show()