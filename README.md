# Tracking Code 
Tracking code for extreme precipitation events - contiguous regions of extreme precipitation - on a cube-sphere grid, including demo


This code package identifies and tracks contiguous regions of precipitation above a threshold. 
Starting from a precipitation field at a number of timesteps on a cubed sphere grid (sample data included, `cubesphere305small.mat`), 
a threshold for precipitation deemed "extreme" is calculated over the entire dataset (`calc95thpercentile.m`), 
a connectivity matrix is for the [NCAR's CESM](http://www2.cesm.ucar.edu/) cube-sphere grid is calculated (`cubesphereneighbors.m`),
applied to the field to make it binary, contiguous regions above the threshold are identified at each timestep (`spatialregionneighborsearchcentroid.m`), 
and these contiguous regions are tracked from one time to the next (`trackcodematrixempty.m`). 
 
Some scientific considerations and an application of the tracking code are shown in the following manuscript, which
you should cite if you use it the code or algorithm: 
<a href="http://www.atmos.washington.edu/~angie/pendergrassreedmedeiros_submitted.pdf">Pendergrass, A.G., K.A. Reed and B. Medeiros, The link between extreme precipitation and convective organization in a warming climate, Submitted to Geophysical Research Letters on 29 June 2016.</a>
