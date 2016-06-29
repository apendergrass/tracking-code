# Tracking Code Demo
Tracking code for extreme precipitation events - contiguous regions of extreme precipitation - on a cube-sphere grid including demo




This code package identifies and tracks contiguous regions of precipitation above a threshold. 
Starting from a precipitation field at a number of timesteps on a cubed sphere grid (sample data included), 
a threshold for precipitation deemed "extreme" is calculated over the entire dataset,
applied to the field to make it binary, contiguous regions above the threshold are identified at each timestep, 
and these contiguous regions are tracked from one time to the next. 
 
Some scientific considerations and an application of the tracking code are shown in the following manuscript, which
you should cite if you use it the code or algorithm: 
Pendergrass, A.G., K.A. Reed and B. Medeiros, The link between extreme precipitation and convective organization in a warming climate, Submitted to Geophysical Research Letters on 29 June 2016. 
