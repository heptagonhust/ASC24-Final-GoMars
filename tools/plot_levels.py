#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import xarray as xr

parser = argparse.ArgumentParser(description='Plot vertical levels')
parser.add_argument('file', type=str, help='Input file')
parser.add_argument('--lat', type=float, default=-24.5, help='Latitude')
args = parser.parse_args()

f = xr.open_dataset(args.file, decode_times=False)

ax = plt.subplot(111)
ax.set_title(f'Vertical levels at latitude {args.lat}')
ax.set_xlabel('Longitude (deg)')
ax.set_xlim(0, 360)
ax.set_ylabel('Height (m)')
ax.set_ylim(0, 50000)
for k in range(f['ilev'].size):
    ax.plot(f['lon'], f['z'].sel(time=0, lat=args.lat, ilev=f['ilev'][k], method='nearest'), linewidth=0.5, color='k')

plt.show()
