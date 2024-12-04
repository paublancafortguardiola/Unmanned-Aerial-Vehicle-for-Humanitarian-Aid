# https://medium.com/spatial-data-science/digital-elevation-model-dem-in-python-758f0ede3af8
# EU users who use the Copernicus DEM in their research are requested to use the following DOI when citing the data source in their publications:
# https://portal.opentopography.org/raster?opentopoID=OTSDEM.032021.4326.3

import math 
import numpy as np
import rasterio as rio
from rasterio.plot import show
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


zona_a_visualitzar = "mont_perdut2"

dem = rio.open(zona_a_visualitzar + ".tif")

dem_array = dem.read(1).astype('float64')

start_i, start_j = 0, 0
#dem_array = dem_array[start_i:start_i + 100, start_j:start_j + 100]

dem_array = dem_array[::-1,:] # we invert the elements on the first dimension to have the plot pointing to the north

R = 6.37*10**6
theta = np.radians(47)

grid_spacing_i = np.radians(1/3600)
grid_spacing_j = 0

if np.radians(40) < theta < np.radians(90):
    grid_spacing_j = np.radians(1/3600)
elif np.radians(20) < theta < np.radians(40):
    grid_spacing_j = np.radians(2/3600)
elif np.radians(15) < theta < np.radians(20):
    grid_spacing_j = np.radians(3/3600)
elif np.radians(10) < theta < np.radians(15):
    grid_spacing_j = np.radians(4/3600)
else:
    grid_spacing_j = np.radians(6/3600)


d_nni = 2*R*np.sin(grid_spacing_i/2)
d_nnj = 2*R*np.sin(theta)*np.sin(grid_spacing_j/2)

 # spacing between points:
# https://spacedata.copernicus.eu/documents/20123/121239/GEO1988-CopernicusDEM-SPE-002_ProductHandbook_I4.0.pdf

# plot a map of the area in which we will work2450
x = np.linspace(0, dem_array.shape[1]*0.001*d_nnj, dem_array.shape[1])
y = np.linspace(0, dem_array.shape[0]*0.001*d_nni, dem_array.shape[0])
X, Y = np.meshgrid(x, y)
Z = dem_array
plt.imshow(Z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
plt.title('Colormap of the altitude')
plt.colorbar()
plt.show()

pendent_permes = np.tan(26/360*2*np.pi)
print(pendent_permes)
print(dem_array.shape)

altitud_maxima = np.max(dem_array)
print(altitud_maxima)
altitud_minima = np.min(dem_array)
print(altitud_minima)

plt.subplot(2, 1, 1)
# Create a grid of values
x = np.linspace(0, dem_array.shape[1]*0.001*d_nnj, dem_array.shape[1])
y = np.linspace(0, dem_array.shape[0]*0.001*d_nni, dem_array.shape[0])
X, Y = np.meshgrid(x, y)
Z = dem_array
plt.imshow(Z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
plt.title('Colormap of the altitude')
plt.colorbar()

# Cone matrix
Dk_max = math.trunc((altitud_maxima-altitud_minima)/d_nni/pendent_permes) + 1
Dl_max = math.trunc((altitud_maxima-altitud_minima)/d_nnj/pendent_permes) + 1

cone_max = np.zeros((2*Dk_max, 2*Dl_max))
for k in range(0, 2*Dk_max):
    for l in range(0, 2*Dl_max):
        cone_max[k,l] = (((k - Dk_max)**2*d_nni**2+(l - Dl_max)**2*d_nnj**2)**0.5)*pendent_permes

dem_array.size
altures_minimes_de_vol = np.zeros(dem_array.shape)
for j in range(0, dem_array.shape[1]):
    for i in range(0, dem_array.shape[0]):

        altitudij = dem_array[i][j]

        # dimensions of the searching zone
        Dk = math.trunc((altitud_maxima-altitudij)/d_nni/pendent_permes) + 1
        Dl = math.trunc((altitud_maxima-altitudij)/d_nnj/pendent_permes) + 1

        k_min = max(0, (i-Dk))
        k_max = min(dem_array.shape[0], (i+Dk))
        l_min = max(0, (j-Dl))
        l_max = min(dem_array.shape[1],j+Dl)

        if k_max-k_min + l_max - l_min > 0:
            altures_minimes_de_vol[i][j] = np.max(dem_array[k_min: k_max, l_min : l_max] - cone_max[(Dk_max -(i - k_min)): (Dk_max+(k_max - i)), (Dl_max-(j - l_min)) : (Dl_max + (l_max - j))] - np.full((k_max-k_min, l_max - l_min), altitudij))
    print(j)

    if j in np.round(np.linspace(0, dem_array.shape[1] -1, 10)).astype(int):
        plt.subplot(2, 1, 2)
        Z = altures_minimes_de_vol
        plt.imshow(Z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
        plt.title('Colormap of the minimum flying altitude over the ground')
        plt.colorbar()
        plt.show(block = False)
        plt.pause(1)

np.save(zona_a_visualitzar + ".npy", altures_minimes_de_vol)
import numpy as np

# % of landing surface
print(f'The percentage of elements with a value of 0 is {(np.count_nonzero(altures_minimes_de_vol.flatten() < 5) / altures_minimes_de_vol.flatten().size) * 100:.2f}%')

plt.show()