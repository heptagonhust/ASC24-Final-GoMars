#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import xarray as xr

parser = argparse.ArgumentParser(description='Plot vertical levels')
parser.add_argument('file', type=str, help='Input file')
parser.add_argument('--lat', type=float, default=-24.5, help='Latitude')
parser.add_argument('--lon1', type=float, default=0, help='Minimum longitude')
parser.add_argument('--lon2', type=float, default=360, help='Maximum longitude')
args = parser.parse_args()

f = xr.open_dataset(args.file, decode_times=False)

lon = f.lon.sel(lon=slice(args.lon1, args.lon2))
z = f.z.sel(time=0, lat=args.lat, lon=lon, method='nearest')

ax = plt.subplot(111)
ax.set_title(f'Vertical levels at latitude {args.lat}')
ax.set_xlabel('Longitude (deg)')
ax.set_xlim(args.lon1, args.lon2)
ax.set_ylabel('Height (m)')
ax.set_ylim(0, z.max())

if 'lev' in f.z.dims:
	for k in range(f.lev.size):
		ax.plot(lon, z.sel(lev=f['lev'][k], method='nearest'), linewidth=0.5, color='k')
else:
	for k in range(f['ilev'].size):
		ax.plot(lon, z.sel(ilev=f['ilev'][k], method='nearest'), linewidth=0.5, color='k')

plt.show()
