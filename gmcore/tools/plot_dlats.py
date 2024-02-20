#!/usr/bin/env python3

import argparse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

parser = argparse.ArgumentParser(description='Set meridional intervals (dlat) and plot them.')
parser.add_argument('--nc')
parser.add_argument('--nlat', type=int)
parser.add_argument('--dlat', type=float)
parser.add_argument('--pole-mul', type=float)
parser.add_argument('--pole-decay', type=float)
args = parser.parse_args()

if args.nc:
	f = Dataset(args.nc, 'r')
	full_nlat = f.variables['lat'].size
	half_nlat = f.variables['ilat'].size
	full_lat = np.radians(f.variables['lat'])
	half_lat = np.radians(f.variables['ilat'])
	dlat = np.diff(full_lat)
else:
	full_nlat = args.nlat
	half_nlat = full_nlat - 1
	dlat0        = np.radians(args.dlat)
	pole_mul     = args.pole_mul
	pole_decay   = args.pole_decay

	full_lat = np.ndarray([full_nlat])
	half_lat = np.ndarray([half_nlat])
	dlat = np.ndarray([half_nlat])

	dlat1 = np.pi / half_nlat
	for j in range(half_nlat):
		half_lat[j] = -np.pi / 2 + (j + 0.5) * dlat1

	for j in range(half_nlat):
		dlat[j] = dlat0 * (1 + (pole_mul - 1) * np.exp(-pole_decay * (np.abs(half_lat[j]) - np.pi / 2)**2))
	dlat = dlat / np.sum(dlat) * np.pi

	half_lat[0] = -np.pi / 2 + dlat[0] / 2
	for j in range(1, half_nlat):
		half_lat[j] = half_lat[j-1] + (dlat[j-1] + dlat[j]) / 2

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_title(f'NLAT = {full_nlat}')
plt.plot(np.degrees(half_lat), np.degrees(dlat))
if 'dlat0' in locals():
	dlat[:] = dlat0
	plt.plot(np.degrees(half_lat), np.degrees(dlat))
plt.show()
plt.close()
