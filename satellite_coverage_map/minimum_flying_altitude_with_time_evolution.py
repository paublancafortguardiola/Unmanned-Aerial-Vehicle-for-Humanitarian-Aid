import numpy as np
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
import math
import rasterio as rio
from rasterio.plot import show
np.set_printoptions(threshold=np.inf)
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter
import multiprocessing
import os
import time

def cartesian_position(radius, theta, phi):
    # function that returns the cartesian position from the angular coordinates of a point
    return radius * np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])

def calculate_cell_size(theta_A):
    # function that calculated the size of the cells: the distance to the nearest neighbor
    grid_spacing_i = np.radians(1/3600)
    grid_spacing_j = 0

    if np.radians(40) < theta_A < np.radians(90):
        grid_spacing_j = np.radians(1/3600)
    elif np.radians(20) < theta_A < np.radians(40):
        grid_spacing_j = np.radians(2/3600)
    elif np.radians(15) < theta_A < np.radians(20):
        grid_spacing_j = np.radians(3/3600)
    elif np.radians(10) < theta_A < np.radians(15):
        grid_spacing_j = np.radians(4/3600)
    else:
        grid_spacing_j = np.radians(6/3600)

    d_nni = 2*R*np.sin(grid_spacing_i/2)
    d_nnj = 2*R*np.sin(theta_A)*np.sin(grid_spacing_j/2)
    return[d_nni,d_nnj]

def satellite_position_in_A_reference_frame(theta_A, phi_A, theta_S_1, phi_S_1):
    # position of the satellite using a local base, with the directions {up, east, south}   
    
    # cartesian position of points A and S
    A = cartesian_position(R, theta_A, phi_A)
    S = cartesian_position(R+h, theta_S_1, phi_S_1)

    # vector going from A to S
    AS = S-A

    # A local base, with the directions {up, east, south}
    rotation_matrix_1 = np.array([[np.cos(phi_A), np.sin(phi_A), 0],[-np.sin(phi_A), np.cos(phi_A), 0], [0,0,1]])
    rotation_matrix_2 = np.array([[np.cos(theta_A), 0, -np.sin(theta_A)], [0,1,0] ,[np.sin(theta_A), 0, np.cos(theta_A)]])
    rotation_matrix = rotation_matrix_2@rotation_matrix_1
    AS = rotation_matrix@AS
    return AS

def calculate_elevation(AS):
    # calculation of the elevation angle that forms AS vector with the xy plane
    elevation = np.arcsin(AS[2]/np.linalg.norm(AS))
    return elevation

def resize_matrix(matrix):
    # resize the matrix
    start_i, start_j = 0, 0
    matrix = matrix[start_i:start_i + 20, start_j:start_j + 20]
    return matrix
    
def preset_dem(zone):
    # funtion that prepares the dem map for its use and find necessaryy properties
    dem = rio.open(zone + ".tif")
    dem_array = dem.read(1).astype('float64')
    # dem_array = resize_matrix(dem_array) # if we are only interested in one subregion, we can select it here
    return dem_array

def dem_infleunce_parameters(dem_array, elevation):
    # calculation of the maximum altitude and the influence distance. The later is defined as the distance at which a mountaing that obstructs the A-Satellite line of sight
    maximum_altitude = np.max(dem_array)
    minimum_altitude = np.min(dem_array)
    influence_distance = (maximum_altitude-minimum_altitude)/np.tan(elevation)
    return [maximum_altitude, influence_distance]

def construct_auxiliar_matrices(Dk_max,Dl_max,elevation,AS,d_nni,d_nnj):
    # function that builds ones_and_zeros_max and cone_max, two necessary matrices to compute the minimum flying altitude matrix
    ones_and_zeros_max = np.zeros((2*Dk_max+1, 2*Dl_max+1))
    cone_max = np.zeros((2*Dk_max+1, 2*Dl_max+1))
    for k in range(2*Dk_max+1):
        if np.sign(k-Dk_max) == np.sign(AS[0]):
            l = Dl_max + round((k-Dk_max)*d_nni*AS[1]/(AS[0]*d_nnj))
            if l in range(2*Dl_max+1):
                ones_and_zeros_max[k,l] = 1 
                cone_max[(k,l)] = (((k-Dk_max)**2*d_nni**2+(l-Dl_max)**2*d_nnj**2)**0.5)*np.tan(elevation)
    for l in range(2*Dl_max+1):
        if np.sign(l-Dl_max) == np.sign(AS[1]):
            k = Dk_max + round((l-Dl_max)*d_nnj*AS[0]/(AS[1]*d_nni))
            if k in range(2*Dk_max+1):
                ones_and_zeros_max[k,l] = 1 
                cone_max[(k,l)] = (((k-Dk_max)**2*d_nni**2+(l-Dl_max)**2*d_nnj**2)**0.5)*np.tan(elevation)
    return [ones_and_zeros_max,cone_max]

def calculate_minimum_flying_altitude(dem_array, ones_and_zeros_max,cone_max, maximum_altitude, elevation, Dk_max, Dl_max,d_nni,d_nnj):
    # calculation of the minimum flying altitude
    minimum_flying_altitude = np.zeros(dem_array.shape)
    for i in range(0, dem_array.shape[0]):
        for j in range(0, dem_array.shape[1]):
            altitudeij = dem_array[i][j]

            # dimensions of the searching zone
            Dk = math.trunc((maximum_altitude-altitudeij)/d_nni/np.tan(elevation)) + 1
            Dl = math.trunc((maximum_altitude-altitudeij)/d_nnj/np.tan(elevation)) + 1

            k_min = max(0, (i-Dk))
            k_max = min(dem_array.shape[0]-1, i+Dk)
            l_min = max(0, (j-Dl))
            l_max = min(dem_array.shape[1]-1,j+Dl)

            if (k_max-k_min) > 0 and (l_max-l_min) > 0:
                minimum_flying_altitude[i][j] = np.max((dem_array[k_min : k_max + 1, l_min : l_max +1] - cone_max[(Dk_max -(i - k_min)): (Dk_max+(k_max - i)) + 1, (Dl_max-(j - l_min)) : (Dl_max + (l_max - j)) + 1] - np.full((k_max-k_min +1, l_max - l_min +1), altitudeij))*ones_and_zeros_max[(Dk_max -(i - k_min)): (Dk_max+(k_max - i) + 1), (Dl_max-(j - l_min)) : (Dl_max + (l_max - j)) + 1])
            #print(i)
    return minimum_flying_altitude


def plot(dem_array, minimum_flying_altitude, d_nni, d_nnj):
    # plot the altutude of the terrain (dem_array) and the minimum flying altitude
    # we reorient the maps
    dem_array = dem_array[::-1,:] 
    minimum_flying_altitude = minimum_flying_altitude[::-1,:] 

    plt.subplot(2, 1, 1)
    x = np.linspace(0, dem_array.shape[1]*d_nnj, dem_array.shape[1])
    y = np.linspace(0, dem_array.shape[0]*d_nni, dem_array.shape[0])
    X, Y = np.meshgrid(x, y)
    Z = dem_array
    plt.imshow(Z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar()

    plt.subplot(2, 1, 2)
    Z = minimum_flying_altitude
    plt.imshow(Z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
    plt.colorbar()
    plt.show()

def plot_video(dem_array, minimum_flying_altitude, d_nni, d_nnj, times, number_of_instants_calculated,number_of_instants_interpolated, output_file='animation.mp4'):
    # we reorient the maps
    dem_array = dem_array[::-1,:] 
    for i, e in enumerate(minimum_flying_altitude):
        minimum_flying_altitude[i] = e[::-1,:] 

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Initial setup
    x = np.linspace(0, dem_array.shape[1] * d_nnj, dem_array.shape[1])
    y = np.linspace(0, dem_array.shape[0] * d_nni, dem_array.shape[0])
    X, Y = np.meshgrid(x, y)

    # DEM plot
    Z1 = dem_array
    im1 = ax1.imshow(Z1, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
    fig.colorbar(im1, ax=ax1)
    ax1.set_title('Altitude above mean sea level')

    # Minimum Flying Altitude plot
    Z2 = minimum_flying_altitude[0]
    im2 = ax2.imshow(Z2, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='viridis', aspect='auto')
    fig.colorbar(im2, ax=ax2)
    ax2.set_title('Minimum Flying Altitude')

    def update(frame):
        Z2 = minimum_flying_altitude[frame]
        im2.set_data(Z2)
        ax2.set_title(f'Minimum Flying Altitude (time: {int(times[0]+frame/(number_of_instants_calculated*number_of_instants_interpolated)*(max(times)-min(times)))}s, Landing area: {int(calculate_landing_area_percentage(minimum_flying_altitude[frame]))}% ).')
        print("frame:")
        print(frame)
        return im2,

    anim = FuncAnimation(fig, update, frames=minimum_flying_altitude.shape[0], blit=True)
    
    anim.save(output_file, writer='ffmpeg')
    print(f"Animation saved as {output_file}")

def plot_3d_video(dem_array, minimum_flying_altitude, d_nni, d_nnj, interval_of_time_calculated,number_of_instants_calculated,number_of_instants_interpolated, output_file='3d_animation.mp4'):
    # we reorient the maps
    dem_array = dem_array[::-1,:] 
    for i, e in enumerate(minimum_flying_altitude):
        minimum_flying_altitude[i] = e[::-1,:] 
 
    minimum_flying_altitude_normalized = (minimum_flying_altitude - np.min(minimum_flying_altitude)) / (np.max(minimum_flying_altitude) - np.min(minimum_flying_altitude))

    fig = plt.figure(figsize=(20, 16))
    ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(0, dem_array.shape[1] * d_nnj, dem_array.shape[1])
    y = np.linspace(0, dem_array.shape[0] * d_nni, dem_array.shape[0])
    X, Y = np.meshgrid(x, y)
    cmap = plt.cm.viridis_r

    norm = plt.Normalize(vmin=np.min(minimum_flying_altitude), vmax=np.max(minimum_flying_altitude))
    mappable = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array([]) 

    def update(frame):
        ax.clear()
        ax.plot_surface(X, Y, dem_array, facecolors=cmap(minimum_flying_altitude_normalized[frame]), rstride=1, cstride=1, linewidth=0, antialiased=False)
        
        ax.set_xlabel('X (meters)')
        ax.set_ylabel('Y (meters)')
        ax.set_zlabel('Altitude (meters)')
        ax.set_zlim(np.max(dem_array), min(d_nni * dem_array.shape[0], d_nnj * dem_array.shape[1]))
        ax.set_title(f'time: {times[0]+frame/(number_of_instants_calculated*number_of_instants_interpolated)*(max(times)-min(times)):.2f}s   Landing area: {calculate_landing_area_percentage(minimum_flying_altitude[frame])}%', y=0.05)
        ax.set_axis_off()
        mappable.set_array(minimum_flying_altitude[frame])
        print("frame:")
        print(frame)
        ax.view_init(elev=30, azim=-90+frame * 0.4)

    anim = FuncAnimation(fig, update, frames=minimum_flying_altitude.shape[0], blit=False)

    writer = FFMpegWriter(fps=10, metadata=dict(artist='Me'), bitrate=1800)

    anim.save(output_file, writer=writer)
    print(f"Animation saved as {output_file}")

def main(theta_A,phi_A,theta_S_1,phi_S_1,dem_array):

    # given a position of the observer, a position of a satellite, and a terrain altitude (dem_array) it returns the minimum flying altitude of the UAV that assures a nonobstructed line of sight. If the Earth's curvature is obstructing the line of sight, the function returns a minimum flying altitude of 10000

    AS = satellite_position_in_A_reference_frame(theta_A, phi_A, theta_S_1, phi_S_1)
    elevation = calculate_elevation(AS)

    if elevation > np.radians(5):
        #print('elevation angle = ', np.degrees(elevation),"º",)
        [maximum_altitude, influence_distance] = dem_infleunce_parameters(dem_array, elevation)

        # calculation of the size of the affecting zone
        Dk_max = math.trunc(influence_distance/d_nni) + 1
        Dl_max = math.trunc(influence_distance/d_nnj) + 1

        [ones_and_zeros_max,cone_max] = construct_auxiliar_matrices(Dk_max,Dl_max,elevation,AS,d_nni,d_nnj)

        minimum_flying_altitude = calculate_minimum_flying_altitude(dem_array, ones_and_zeros_max,cone_max, maximum_altitude, elevation, Dk_max, Dl_max, d_nni,d_nnj)
        return minimum_flying_altitude
    else:
        return 10000*np.ones(dem_array.shape)

def linear_interpolation (A_matrix,B_matrix, n_interpolation_points):
    # function that returns n_interpolation_points matrices corresponding to the linear interpolation between matrix A and matrix B.
    parameters = np.linspace(0, 1-1/n_interpolation_points, n_interpolation_points)
    interpolated_matrix = np.array([np.zeros((A_matrix.shape[0], A_matrix.shape[1])) for _ in range(n_interpolation_points)])
    for index, parameter in enumerate(parameters):
        interpolated_matrix[index] = (1-parameter)*A_matrix + parameter*B_matrix
    return interpolated_matrix

def calculate_minimum_flying_altitude_for_a_window(dem_array, theta_A, phi_A, satellites, zone, times, CPUs_used = 1):
    # function that calculates and plots the minimum flying altitude in a given terrain (dem_arrauy), situated at a certain position (theta_A, phi_A,), connected to a constellation of satellites (satellites).
    
    # angular speed of the satellites and the observer A. These values are valid for the Iridium constellation
    omega_S = - 2*np.pi/(100*60)
    omega_A = 2*np.pi/(24*60*60)
    
    phi_A = phi_A*np.ones(times.shape) + times*omega_A
    
    for i, satellite in enumerate(satellites):
        satellite.set_theta(satellite.theta_initial*np.ones(times.shape) + times*omega_S)

    matrices_index = range(len(times))
    minimum_flying_altitude_calculated = np.array([np.ones((dem_array.shape[0], dem_array.shape[1])) for _ in matrices_index])

    args = [(instant, dem_array, theta_A, phi_A, satellites) for instant in matrices_index]

    if CPUs_used == 1:
        for instant in matrices_index:
            minimum_flying_altitude_calculated[instant] = minimum_flying_altitude_for_all_satellites(instant, dem_array, theta_A, phi_A, satellites)
    else: 
        with multiprocessing.Pool(CPUs_used) as pool: 
            minimum_flying_altitude_calculated = pool.starmap(minimum_flying_altitude_for_all_satellites, args)
    return minimum_flying_altitude_calculated

def minimum_flying_altitude_for_all_satellites(instant, dem_array, theta_A, phi_A, satellites):
    start_time = time.time()
    minimum_flying_altitude = 10000*np.ones(dem_array.shape)
    print("time of the calculation = ", times[instant],"s started")
    for i, satellite in enumerate(satellites):
        #print("")
        #print("time = ", times[instant],"s")
        minimum_flying_altitude_Si = main(theta_A,phi_A[instant],satellite.theta[instant],satellite.phi_initial,dem_array)
        minimum_flying_altitude = np.minimum(minimum_flying_altitude_Si,minimum_flying_altitude)
    print()
    print("time = ", times[instant],"s finished")
    print(f"Execution time: {time.time()-start_time:.4f} seconds")
    print()
    return minimum_flying_altitude

def calculate_landing_area_percentage(minimum_flying_altitude):
    # function that calculates the landing area in a minimum_flying_altitude We define the landing altitude as the area in which the minimum flying altitude is less than 5 meters.
    return  (np.count_nonzero(minimum_flying_altitude.flatten() < height_threshold) / minimum_flying_altitude.flatten().size) * 100

class Satellite:
    # class that will store the satellite data
    def __init__(self, phi_initial, theta_initial):
        self.phi_initial = phi_initial
        self.theta_initial = theta_initial
        self.theta = None

    def set_theta(self, theta):
        self.theta = theta

def satellites_angular_position():
    # angular position of all the satellites of the Iridium constellation if we assume that there is a satellite at north pole, and the others are consequently positioned.

    # angular distance between satellites
    Dphi = 31.6 * 2 * np.pi / 360

    satellites = []
    for i in range(66):
        satellites.append(Satellite(Dphi*np.trunc(i/11), np.radians(360)/11*((i)%11) + ((np.trunc(i/11)%11)%2)*np.radians(360)/22))
    return satellites

def interpolate_minimum_flying_altitude_for_a_window(minimum_flying_altitude_calculated,number_of_instants_interpolated,times):
    minimum_flying_altitude_interpolated = np.array([np.zeros((dem_array.shape[0], dem_array.shape[1])) for _ in range((len(times)-1)*number_of_instants_interpolated+1)])
    for instant in range(len(times)-1):
        minimum_flying_altitude_interpolated[instant*number_of_instants_interpolated : (instant+1)*number_of_instants_interpolated] = linear_interpolation(minimum_flying_altitude_calculated[instant], minimum_flying_altitude_calculated[instant+1], number_of_instants_interpolated)
    minimum_flying_altitude_interpolated[-1] = minimum_flying_altitude_calculated[-1]
    return minimum_flying_altitude_interpolated

def plot_landing_surface(minimum_flying_altitude, times):
    landing_surfaces = []
    for i in range(minimum_flying_altitude.shape[0]):
        landing_surfaces.append(calculate_landing_area_percentage(minimum_flying_altitude[i]))

    plt.figure(figsize=(10, 6))
    plt.plot(np.linspace(min(times), max(times), minimum_flying_altitude.shape[0]), landing_surfaces, label='Landing Surface')
    plt.xlabel('Time')
    plt.ylabel('Landing Surface (%)')
    plt.title('Landing Surface as a Function of Time')
    plt.legend()
    plt.grid(True)
    plt.savefig('landing_surface.png') 
    plt.close()

## Programm
# expected computation time of this example: 2h. 
# For quicker answers, chose an other zone, or select a subregion in the function "preset_dem()""

# Constants
R = 6.37e6  # Radius of the Earth in meters
h = 7.8e5   # Altitude of the satellites in meters

# angular position of the observer A. In this case it is the location of the Ordesa and Monte Perdido National Park
theta_A = np.radians(49)
phi_A = np.radians(2) 

# height threshole: if the minimum flying altitude is below this one we will considere that the UAV can land
height_threshold = 5

# terrain imput:
zone = "mont_perdut2"

# preparation for the calculation

# cell size
[d_nni,d_nnj] = calculate_cell_size(theta_A)

# contruction of the terrain altitude of the zone "zone"
dem_array= preset_dem(zone)

# definition of the position of all the satellites constellation
satellites = satellites_angular_position()

print('Number of CPUs in the system: {}'.format(os.cpu_count()))
CPUs_used=1 # to select the one with better performance

time_begin = 0 # time at which the calculation begins in minutes
time_end = 10 # time at which the calculation ends in minutes
number_of_instants_calculated = 11

times = np.linspace(time_begin, time_end, number_of_instants_calculated)*60
minimum_flying_altitude_calculated = calculate_minimum_flying_altitude_for_a_window(dem_array, theta_A, phi_A, satellites, zone,times,CPUs_used)
np.save(zone + ".npy", minimum_flying_altitude_calculated)

number_of_instants_interpolated = 10 # between calculated instants
minimum_flying_altitude_interpolated = interpolate_minimum_flying_altitude_for_a_window(minimum_flying_altitude_calculated,number_of_instants_interpolated,times)

# plot the results
plot_landing_surface(minimum_flying_altitude_interpolated, times)

plot_video(dem_array, minimum_flying_altitude_interpolated, d_nni, d_nnj,times,number_of_instants_calculated,number_of_instants_interpolated)

plot_3d_video(dem_array, minimum_flying_altitude_interpolated, d_nni, d_nnj,times,number_of_instants_calculated,number_of_instants_interpolated)
