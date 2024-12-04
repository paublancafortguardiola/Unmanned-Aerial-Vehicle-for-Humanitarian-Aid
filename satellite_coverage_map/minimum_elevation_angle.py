import numpy as np
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
import math

# Constants
R = 6.37e6  # Radius of the Earth in meters
h = 7.8e5   # Altitude in meters
phi = np.radians(2)     # Initial azimuthal angle
theta = np.radians(47)  # Initial polar angle in radians

def distancia_minima_1_satellits(x):
    phi_rel, theta_rel = x
    A = R * np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    S = (R + h) * np.array([np.cos(phi + phi_rel) * np.sin(theta + theta_rel), 
                            np.sin(phi + phi_rel) * np.sin(theta + theta_rel), 
                            np.cos(theta + theta_rel)])
    return np.linalg.norm(S - A), [A, S]

def distancia_minima_3_satellits(x):
    phi_rel, theta_rel = x
    A = R * np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    P = (R + h) * np.array([np.cos(phi + phi_rel - Dphi) * np.sin(theta + theta_rel + Dtheta / 2), 
                            np.sin(phi + phi_rel - Dphi) * np.sin(theta + theta_rel + Dtheta / 2), 
                            np.cos(theta + theta_rel + Dtheta / 2)])
    Q = (R + h) * np.array([np.cos(phi + phi_rel - Dphi) * np.sin(theta + theta_rel - Dtheta / 2), 
                            np.sin(phi + phi_rel - Dphi) * np.sin(theta + theta_rel - Dtheta / 2), 
                            np.cos(theta + theta_rel - Dtheta / 2)])
    S = (R + h) * np.array([np.cos(phi + phi_rel) * np.sin(theta + theta_rel), 
                            np.sin(phi + phi_rel) * np.sin(theta + theta_rel), 
                            np.cos(theta + theta_rel)])
    return np.min([np.linalg.norm(P - A), np.linalg.norm(Q - A), np.linalg.norm(S - A)]), [A, P, Q, S]

def distancia_minima_4_satellits(x):
    phi_rel, theta_rel1 = x
    theta_rel120 = Dtheta - 2*((theta + Dtheta/4) % (Dtheta/2))

    if theta_rel1 < theta_rel120:
        theta_rel2 = theta_rel120 - theta_rel1
    else:
        theta_rel2 = Dtheta + theta_rel120 - theta_rel1

    A = R * np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    P = (R + h) * np.array([np.cos(phi + phi_rel - Dphi) * np.sin(theta + theta_rel2), 
                            np.sin(phi + phi_rel - Dphi) * np.sin(theta + theta_rel2), 
                            np.cos(theta + theta_rel2)])
    Q = (R + h) * np.array([np.cos(phi + phi_rel - Dphi) * np.sin(theta + theta_rel2 - Dtheta), 
                            np.sin(phi + phi_rel - Dphi) * np.sin(theta + theta_rel2 - Dtheta), 
                            np.cos(theta + theta_rel2 - Dtheta)])
    S = (R + h) * np.array([np.cos(phi + phi_rel) * np.sin(theta + theta_rel1), 
                            np.sin(phi + phi_rel) * np.sin(theta + theta_rel1), 
                            np.cos(theta + theta_rel1)])
    T = (R + h) * np.array([np.cos(phi + phi_rel) * np.sin(theta + theta_rel1 - Dtheta), 
                            np.sin(phi + phi_rel) * np.sin(theta + theta_rel1 - Dtheta), 
                            np.cos(theta + theta_rel1 - Dtheta)])
    return np.min([np.linalg.norm(P - A), np.linalg.norm(Q - A), np.linalg.norm(S - A), np.linalg.norm(T - A)]), [A, P, Q, S, T]

def neg_distancia_minima_3_satellits_boundary(x):
    phi_rel = x[0]
    theta_rel = x[1] * phi_rel / Dphi
    return -distancia_minima_3_satellits([phi_rel, theta_rel])[0]

def neg_distancia_minima_4_satellits_boundary(x):
    return -distancia_minima_4_satellits(x)[0]

def calculate_elevation_angle(result):
    return np.degrees(np.arccos(-((R + h)**2 - R**2 - result.fun**2) / (2 * R * (-result.fun))) - np.pi / 2)

def calculate_elevation_angle2(dist):
    return np.degrees(np.arccos(-((R + h)**2 - R**2 - dist**2) / (2 * R * (dist))) - np.pi / 2)

# Define the function to calculate elevation angle given theta for each case
def calculate_min_elevation_angle(theta, case=1):
    global Dtheta, Dphi
    if case == 1:
        Dtheta = 33 * 2 * np.pi / 360
        Dphi = 31.6 * 2 * np.pi / 360
        bounds = [(0, Dphi), (-Dtheta / 2, Dtheta / 2)]
        result = differential_evolution(neg_distancia_minima_3_satellits_boundary, bounds, tol=1e-8)
    elif case == 2:
        Dtheta = 33 * 2 * np.pi / 360
        Dphi = 22 * 2 * np.pi / 360
        bounds = [(0, Dphi), (0, Dtheta)]
        result = differential_evolution(neg_distancia_minima_4_satellits_boundary, bounds, tol=1e-8)
    
    # Calculate the elevation angle
    elevation = np.arccos(((R + h)**2 - R**2 - (result.fun)**2) / (2 * R * result.fun)) - np.pi / 2
    return elevation * 360 / (2 * np.pi)

# Plot Earth surface as a sphere
def plot_earth_surface(ax):
    definition = 50
    u = np.linspace(-np.pi, np.pi, 2*definition)
    v = np.linspace(0, np.pi, definition)
    x = (R) * np.outer(np.cos(u), np.sin(v))
    y = (R) * np.outer(np.sin(u), np.sin(v))
    z = (R) * np.outer(np.ones(np.size(u)), np.cos(v))

    #Load Earth texture image
    earth_texture = Image.open('earth_texture2.jpg')
    earth_texture = np.array(earth_texture.resize((2*definition, definition))) / 255.0
    earth_texture = np.transpose(earth_texture, (1, 0, 2))

    ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=earth_texture, edgecolor='none', alpha=1, zorder=1)

    ax.view_init(elev=theta*180/np.pi, azim=0)

# Case 1: Three satellites
Dtheta, Dphi = np.radians(33), np.radians(31.6)
bounds = [(0, Dphi), (-Dtheta / 2, Dtheta / 2)]
result = differential_evolution(neg_distancia_minima_3_satellits_boundary, bounds, tol=1e-8)
elevation1 = calculate_elevation_angle(result)
print("Elevation angle for 3 satellites:", elevation1)

# Retrieve positions for plotting
_, satellite_positions_3 = distancia_minima_3_satellits(result.x)

# Case 2: Four satellites
Dphi = np.radians(22)
bounds = [(0, Dphi), (0, Dtheta)]
result = differential_evolution(neg_distancia_minima_4_satellits_boundary, bounds, tol=1e-8)
elevation2 = calculate_elevation_angle(result)
print("Elevation angle for 4 satellites:", elevation2)

# Retrieve positions for plotting
_, satellite_positions_4 = distancia_minima_4_satellits(result.x)

# Final minimum elevation angle
elevation = min(elevation1, elevation2)
print("Minimum elevation angle:", elevation)

# Plotting the 3-satellite case
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d', computed_zorder=False)
ax1.set_title("3-Satellite Configuration")
plot_earth_surface(ax1)

# Plot observer position (A)
A_pos_3 = satellite_positions_3[0]
ax1.scatter(*A_pos_3, color='cyan', s=10, zorder=10, label="Observer A")

# Plot satellites and add to legend
ax1.scatter(*satellite_positions_3[1], color='red', s=10, zorder=10, label="Satellite P")  # P
ax1.scatter(*satellite_positions_3[2], color='orange', s=10, zorder=10, label="Satellite Q")  # Q
ax1.scatter(*satellite_positions_3[3], color='yellow', s=10, zorder=10, label="Satellite S")  # S

# Plot settings
ax1.set_xlabel("X (meters)")
ax1.set_ylabel("Y (meters)")
ax1.set_zlabel("Z (meters)")
ax1.legend()

# Plotting the 4-satellite case
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d', computed_zorder=False)
ax2.set_title("4-Satellite Configuration")
plot_earth_surface(ax2)

# Plot observer position (A)
A_pos_4 = satellite_positions_4[0]
ax2.scatter(*A_pos_4, color='cyan', s=10, zorder=10, label="Observer A")

# Plot satellites and add to legend
ax2.scatter(*satellite_positions_4[1], color='red', s=10, zorder=10, label="Satellite P")  # P
ax2.scatter(*satellite_positions_4[2], color='orange', s=10, zorder=10, label="Satellite Q")  # Q
ax2.scatter(*satellite_positions_4[3], color='yellow', s=10, zorder=10, label="Satellite S")  # S
ax2.scatter(*satellite_positions_4[4], color='green', s=10, zorder=10, label="Satellite T")  # T

# Plot settings
ax2.set_xlabel("X (meters)")
ax2.set_ylabel("Y (meters)")
ax2.set_zlabel("Z (meters)")
ax2.legend()

#plt.show()

# Generate values for theta and corresponding minimum elevation angle
theta_values = np.linspace(2/9*np.pi / 2, np.pi / 2, 2)  # 0 to 90 degrees in radians
min_elevation_angles1 = []
min_elevation_angles2 = []
min_elevation_angles = []

# Calculate the minimum elevation angle between the two cases for each theta
for theta in theta_values:
    elevation_case1 = calculate_min_elevation_angle(theta, case=1)
    elevation_case2 = calculate_min_elevation_angle(theta, case=2)
    min_elevation_angles1.append(elevation_case1)
    min_elevation_angles2.append(elevation_case2)
    min_elevation_angles.append(min(elevation_case1, elevation_case2))

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(theta_values * 180 / np.pi, min_elevation_angles1, color='orange', label='Minimum Elevation Angle for 3-Satellite Configuration')
plt.plot(theta_values * 180 / np.pi, min_elevation_angles2, color='blue', label='Minimum Elevation Angle for 4-Satellite Configuration')
plt.plot(theta_values * 180 / np.pi, min_elevation_angles, color='green', label='Minimum Elevation Angle for both cases')
plt.xlabel('Theta (degrees)')
plt.ylabel('Minimum Elevation Angle (degrees)')
plt.title('Minimum Elevation Angle vs. Theta')
plt.legend()
plt.grid(True)
#plt.show()

N = 501
theta = np.radians(47)
elevation_angles = np.zeros([N,N])
Dtheta, Dphi = np.radians(33), np.radians(31.6)
phi_rel  = np.linspace(0, Dphi, N)
theta_rel = np.linspace(-Dtheta / 2, Dtheta / 2 , N)
for i in range(N):
    for j in range(N):
        if  ((N-1)/2 - i/2)  <= j <= ((N-1)/2 + i/2):         
            dist = distancia_minima_3_satellits([phi_rel[i],theta_rel[j]])[0]+1
            elevation_angles[i][j] = calculate_elevation_angle2(dist)

#elevation_angles = sorted(elevation_angles)
#print(elevation_angles[math.trunc(len(elevation_angles)/2)])

fig4 = plt.figure()
plt.imshow(np.transpose(elevation_angles[::-1,:]), extent=[0, Dphi, -Dtheta / 2, Dtheta / 2], origin='lower', aspect='auto')
plt.colorbar(label='Elevation Angle')
plt.title('Elevation Angles')
plt.xlabel('Phi Relative (radians)')
plt.ylabel('Theta Relative (radians)')
#plt.show()


theta = np.radians(60)
elevation_angles = np.zeros([N,N])
Dtheta, Dphi = np.radians(33), np.radians(31.6)
phi_rel  = np.linspace(-Dphi/2, Dphi/2, N)
theta_rel = np.linspace(-Dtheta / 2, Dtheta / 2 , N)
for i in range(N):
    for j in range(N):        
            dist = distancia_minima_1_satellits([phi_rel[i],theta_rel[j]])[0]+1
            elevation_angles[i][j] = calculate_elevation_angle2(dist)

fig5 = plt.figure()
plt.imshow(np.transpose(elevation_angles[::-1,:]), extent=[-Dphi/2, Dphi/2, -Dtheta / 2, Dtheta / 2], origin='lower', aspect='auto')
plt.colorbar(label='Elevation Angle')
plt.title('Elevation Angles')
plt.xlabel('Phi Relative (radians)')
plt.ylabel('Theta Relative (radians)')
plt.show()

periode = 100.3*60
durada_finestra = 6*60
velocitat_angular = 2*np.pi/periode
mida_finestra = velocitat_angular * durada_finestra
print(mida_finestra)
dist = distancia_minima_1_satellits([0,mida_finestra/2])[0]
elevation_angle = calculate_elevation_angle2(dist)
print(elevation_angle)