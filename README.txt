# CLustre: semi-automated lineament clustering for palaeo-glacial reconstruction
-------------

This script accompanies the paper of Smith et al. (In Press) in ESPL. Here we present a semi-automated algorithm, CLustre, for lineament clustering that uses a locally adaptive, region growing, methodology. CLustre is demonstrated with a polyline data set (ESRI shapefile) representing glacial landforms for palaeo-glacial reconstruction. 

NOTE: This script is accompanied by example data from Stokes & Clark, 2003. When used please cite: 
Stokes, C.R. & Clark, C.D. (2003). The Dubawnt Lake palaeo-ice stream: evidence for dynamic ice sheet behaviour on the Canadian Shield and insights regarding the controls on ice-stream location and vigour. Boreas 32: 263-279. doi: 10.1111/j.1502-3885.2003.tb01442.x

Input: shapefile
Output: edited shapefile with new attribute fields and values
Usage: run Clustre.py with arguments from bash/terminal 

-------------
Arguments:
* lines='location/to/shapefile.shp'
* th_distance=threshold distance
* th_orient=threshold orientation
* th_length=threshold length
* label=classification label
* seed_ID=FID of initial seed lineament
* plotYN='yes' or 'no' to plot results in preview.png

example:

in linux: python Clustre.py lines=\'example/dubawnt.shp\' th_distance=10000 th_orient=2 th_length=3000 label=1 seed_ID=5175 plotYN=\'yes\'
in windows: python Clustre.py lines='example/dubawnt.shp' th_distance=10000 th_orient=2 th_length=3000 label=1 seed_ID=5715 plotYN='yes'

-------------
dependencies (and tested version numbers in brackets): 
* python 2.7        (2.7.9)
* gdal/ogr >= 1.11  (1.11.2)
* numpy             (1.9.2)
* matplotlib        (1.4.3)

Details are described in Smith MJ, Anders NS, Keesstra SD, In Press. CLustre: semi-automated lineament clustering for palaeo-glacial reconstruction. Earth Surface Processes and Landforms

For further questions please contact Niels Anders: n.s.anders@uva.nl

(c) Niels Anders, Mike Smith
Computational Geo-Ecology, Institute for Biodiversity and Ecosystem Dynamics
University of Amsterdam, The Netherlands

