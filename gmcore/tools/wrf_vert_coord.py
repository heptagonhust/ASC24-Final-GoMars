#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

nlev        = 64
t0          = 290           # K
p0          = 100000        # Pa
ptop        = 100           # Pa
dzmax       = 1000          # m
dzbot       = 50            # m
dzstretch_u = 1.1
dzstretch_s = 1.3
rd          = 287           # J/kg/K
g           = 9.81          # m/s2

# Is ztop reasonable with isothermal temperature?
ztop = rd * t0 / g * np.log(p0 / ptop) # m
pbot = p0 * np.exp(-g * dzbot / rd / t0) # Pa

# Allocate arrays.
eta_lev = np.ndarray(nlev+1)
eta = np.ndarray(nlev)
z_lev = np.ndarray(nlev+1)
p_lev = np.ndarray(nlev+1)

# Bottom-up order
eta_lev[0] = 1
eta_lev[1] = (pbot - ptop) / (p0 - ptop)
z_lev[0] = 0
z_lev[1] = dzbot
p_lev[0] = p0
p_lev[1] = pbot

dz = dzbot
for k in range(2, nlev + 1):
	dz = dz * dzstretch_u + (dzstretch_s - dzstretch_u) * np.max([(dzmax / 2 - dz) / (dzmax / 2), 0])
	if (ztop - z_lev[k-1]) / (nlev - k + 1) < dz:
		k -= 1
		break
	z_lev[k] = z_lev[k-1] + dz
	p_lev[k] = p0 * np.exp(-g * z_lev[k] / rd / t0)
	eta_lev[k] = (p_lev[k] - ptop) / (p0 - ptop)
	if k == nlev:
		print('ERROR: Increase nlev!')
		exit(1)
dz = (ztop - z_lev[k]) / (nlev - k)
if dz > 1.6 * dzmax:
	print('ERROR: Upper levels may be too coarse!')
	exit(1)
for k in range(k + 1, nlev + 1):
	z_lev[k] = z_lev[k-1] + dz
	p_lev[k] = p0 * np.exp(-g * z_lev[k] / rd / t0)
	eta_lev[k] = (p_lev[k] - ptop) / (p0 - ptop)
eta_lev[nlev] = 0

dz = np.ndarray(nlev)
dp = np.ndarray(nlev)
deta = np.ndarray(nlev)
for k in range(nlev):
	eta[k] = (eta_lev[k] + eta_lev[k+1]) / 2
	dz[k] = z_lev[k+1] - z_lev[k]
	dp[k] = p_lev[k] - p_lev[k+1]
	deta[k] = eta_lev[k] - eta_lev[k+1]

eta_b = 0.2
eta_c = 0.8
b0 = 2 * eta_b**2 / (1 - eta_b)**3
b1 = -eta_b * (4 + eta_b + eta_b**2) / (1 - eta_b)**3
b2 = 2 * (1 + eta_b + eta_b**2) / (1 - eta_b)**3
b3 = -(1 + eta_b) / (1 - eta_b)**3
c0 = 2 * eta_c**2 / (1 - eta_c)**3
c1 = -eta_c * (4 + eta_c + eta_c**2) / (1 - eta_c)**3
c2 = 2 * (1 + eta_c + eta_c**2) / (1 - eta_c)**3
c3 = -(1 + eta_c) / (1 - eta_c)**3
b = np.zeros(nlev+1)
c = np.zeros(nlev+1)
for k in range(nlev+1):
	if eta_lev[k] >= eta_b:
		b[k] = b0 + b1 * eta_lev[k] + b2 * eta_lev[k]**2 + b3 * eta_lev[k]**3
	if eta_lev[k] >= eta_c:
		c[k] = c0 + c1 * eta_lev[k] + c2 * eta_lev[k]**2 + c3 * eta_lev[k]**3

ps = 65257.574122797363
ref_ps_perb = -1381.6313327777505
p = np.zeros(nlev+1)
for k in range(nlev+1):
	p[k] = b[k] * (ps - ref_ps_perb - ptop) + (eta_lev[k] - b[k]) * (p0 - ptop) + c[k] * ref_ps_perb + ptop

ax = plt.subplot(111)	
ax.plot(p, eta_lev)
ax.invert_yaxis()
ax.set_title('Pressure (Pa)')
plt.show()

# ax1 = plt.subplot(221)
# ax2 = plt.subplot(222)
# ax3 = plt.subplot(223)
# ax4 = plt.subplot(224)

# plt.subplots_adjust(hspace=0.5)
# ax1.plot(z_lev, eta_lev)
# ax1.invert_yaxis()
# ax1.set_title('Height (m)')
# ax2.plot(dz, eta)
# ax2.invert_yaxis()
# ax2.set_title('Height depth (m)')
# ax3.plot(dp, eta)
# ax3.invert_yaxis()
# ax3.set_title('Pressure depth (Pa)')
# ax4.plot(deta, eta)
# ax4.invert_yaxis()
# ax4.set_title('Eta depth (1)')
# plt.show()