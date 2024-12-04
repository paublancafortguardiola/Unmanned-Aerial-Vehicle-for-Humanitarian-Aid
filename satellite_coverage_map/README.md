# UAV Satellite Communication Coverage Map
Satellite communication coverage can be determining for the performance of Beyond Line of Sight (BLOS) Unmanned Aerial Vehicles (UAV). One of the critical problems arises when a mountain obstructs the line of sight between the UAV and the Satellite. 
The goal of this project is to draw a map of the minimum altitude over the ground at which there are no mountains covering the sky above the minimum elevation angle of the satellite.

## Contents
This directory contains the following files:
- "minimum_elevation_angle.py": it calcualtes the minimum elevation angle at which an observer A situated at some given latitude and longitude will find a satellite of the Iridium constellation in the worst-case scenario
- "minimum_flying_altitude.py": given an elevation angle and a terrain DEM map, it calculates a map of the minimmum altitude over the terrain at which the UAV has to fly to be sure that there is clear sky above a the given elevation angle.
- "minimum_flying_altitude_with_time_evolution.py": it is the most complex script. Given a DEM map and an initial position of the satellites of the Iridium constallation calculates the map of the minimum flying altitude at which the UAV has to fly to have a clear line of sight with a satellite, and its time evolution.
- "animation.mp4" and "3d.animation.mp4": videos resulting to use "minimum_flying_altitude_with_time_evolution.py" with "mont_perdut2".

## Example
An example of the Mont Perdut mountain (Pyrenees) is provided, using a DEM map from Copernicus. The next figure shows the resulting 3D map:

![Mont Perdut coverage map](images/mont_perdut.png) 

- **Source**: European Space Agency (2024). *Copernicus Global Digital Elevation Model*. Distributed by OpenTopography. [https://doi.org/10.5069/G9028PQB](https://doi.org/10.5069/G9028PQB). Accessed: 2024-11-03.
