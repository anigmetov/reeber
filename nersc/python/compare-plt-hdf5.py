import numpy as np
import h5py

f = h5py.File('step03387.hdf5', 'r')

# Get constants
universe = f['universe']
hubble = universe.attrs['hubble']
omega_b = universe.attrs['omega_b']
omega_m = universe.attrs['omega_m']
redshift = universe.attrs['redshift']

G = 6.673e-8
rho_c_units_conversion = 1.0e10 * 3.08568025e24 * 5.02785431e-34
pi = np.pi
rho_c = (3e4*hubble**2 / (8 * pi * G)) * rho_c_units_conversion

b = f['/native_fields/baryon_density'][:]
m = f['/native_fields/matter_density'][:]
plt = np.load('plt03387.npy')

bb = b*omega_b*rho_c
mm = m*(omega_m - omega_b)*rho_c

print(np.max(np.abs((bb+mm) - plt)))
