# -*- coding: utf-8 -*-
"""
CLustre: semi-automated lineament clustering for palaeo-glacial reconstruction
-------------
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
in linux:   python Clustre.py lines=\'example_data/example.shp\' th_distance=5000 th_orient=5 th_length=5000 label=1 seed_ID=2511 plotYN=\'yes\'
in windows: python Clustre.py lines='example_data/example.shp' th_distance=5000 th_orient=5 th_length=5000 label=1 seed_ID=2511 plotYN='yes'


-------------
dependencies (and tested version numbers in brackets): 
* python 2.7        (2.7.9)
* gdal/ogr >= 1.11  (1.11.2)
* numpy             (1.9.2)
* matplotlib        (1.4.3)

Details are described in Smith MJ, Anders NS, Keesstra SD, In Press. CLustre: semi-automated lineament clustering for palaeo-glacial reconstruction. Earth Surface Processes and Landforms

For further questions please contact Niels Anders: n.s.anders@uva.nl

(c) Niels Anders
Computational Geo-Ecology, Institute for Biodiversity and Ecosystem Dynamics
University of Amsterdam, The Netherlands
"""

__author__  = 'Niels Anders'
__version__ = '1.0'
__date__    = '2015/08/01'

from osgeo import ogr
import numpy as np
import pylab as plt
import point_store, sys

# ignore warning messages
np.seterr(invalid='ignore')
np.seterr(divide='ignore')

def saveshp(fn, points, proj4):
    point_store.save(fn, points, proj4)
    print 'saved as ',fn   

def readShapefile(lines):
    shapeData = ogr.Open(lines, 1)
    layer = shapeData.GetLayer()
    box = layer.GetExtent()
    nrLines = layer.GetFeatureCount()
    line = layer.GetNextFeature()
       
    print 'nr of features: %2d' % (nrLines)
    # reserve space in memory
    x_start = np.zeros(nrLines)
    y_start = np.zeros(nrLines)
    x_end = np.zeros(nrLines)
    y_end = np.zeros(nrLines)
    length = np.zeros(nrLines)
    x_centroid = np.zeros(nrLines)
    y_centroid = np.zeros(nrLines)
    IDs = np.zeros(nrLines)
    classified = np.zeros(nrLines)
    i = 0
    while line:
        geometry = line.geometry()
        if geometry.GetPointCount() != 2:
            print 'Skipping FID: %5d -> number of points: %2d' % (line.GetFID(), geometry.GetPointCount())
        elif (geometry.GetX(0) == geometry.GetX(1)) & (geometry.GetY(0) == geometry.GetY(1)):
            print 'Skipping FID: %5d -> same end node as start node (zero length)' % line.GetFID()
        else:
            x_start[i]      = geometry.GetX(0)
            x_end[i]        = geometry.GetX(1)
            y_start[i]      = geometry.GetY(0)
            y_end[i]        = geometry.GetY(1)
            x_centroid[i]   = geometry.Centroid().GetX()
            y_centroid[i]   = geometry.Centroid().GetY()
            IDs[i]          = line.GetFID()
            try: classified[i]   = line.GetField('label')
            except: classified[i] = 0
        line = layer.GetNextFeature()
        i += 1
    
    length = np.sqrt((x_end-x_start)**2 + (y_end - y_start)**2)
    ON, OE = CalcOrientation(x_start, x_end, y_start, y_end)
    
    return box, IDs, x_start, x_end, y_start, y_end, x_centroid, y_centroid, length, ON, OE, classified    

def unclassified(IDs, x_start, x_end, y_start, y_end, x_centroid, y_centroid, length, ON, OE, classified, label):
    IDs = IDs[(classified==0) | (classified==label)]
    x_start = x_start[(classified==0) | (classified==label)]
    x_end = x_end[(classified==0) | (classified==label)]
    y_start = y_start[(classified==0) | (classified==label)]
    y_end = y_end[(classified==0) | (classified==label)]
    x_centroid = x_centroid[(classified==0) | (classified==label)]
    y_centroid = y_centroid[(classified==0) | (classified==label)]
    length = length[(classified==0) | (classified==label)]
    ON = ON[(classified==0) | (classified==label)]
    OE = OE[(classified==0) | (classified==label)]
    return IDs, x_start, x_end, y_start, y_end, x_centroid, y_centroid, length, ON, OE
    
    
def CalcOrientation(x_start, x_end, y_start, y_end):
    temp = np.rad2deg(np.arctan((x_end-x_start)/(y_end-y_start)))

    ON = np.zeros(len(temp))
    ON[temp < 0] = 180+(temp[temp<0])
    ON[temp > 0] = temp[temp > 0]

    OE = 180-ON
    return ON, OE

def calcRD(x1, y1, x2, y2):
    dx, dy = x1-x2, y1-y2
    if (dx > 0)&(dy > 0): rd = 90-np.rad2deg(np.arctan(abs(dy/dx)))
    if (dx > 0)&(dy < 0): rd = 90+np.rad2deg(np.arctan(abs(dy/dx)))
    if (dx < 0)&(dy < 0): rd = 270-np.rad2deg(np.arctan(abs(dy/dx)))
    if (dx < 0)&(dy > 0): rd = 270+np.rad2deg(np.arctan(abs(dy/dx)))
    if dx == dy == 0: rd = None
    try:return rd
    except:return None

def labelling(seed_ID, th_length, th_orient, th_distance, fids, FID, ON, OE, X, Y, L, order):
    # seed
    seed_LE, seed_XY, seed_ON, seed_OE = L[seed_ID], (X[seed_ID], Y[seed_ID]), ON[seed_ID], OE[seed_ID]

    # surroundings
    # ... select bounding box
    bb = [seed_XY[0]-th_distance, seed_XY[0]+th_distance, seed_XY[1]-th_distance, seed_XY[1]+th_distance]
    
    # .. search for neighbors
    pneighbors = [(X>bb[0])&(X<bb[1])&(Y>bb[2])&(Y<bb[3])]
    pneighbors = FID[pneighbors]
    for i in pneighbors:            # delete neighbors that have already been used as seed
        if np.where(fids==i)[0].shape[0] > 0: pneighbors = np.delete(pneighbors, pneighbors.tolist().index(i))

    # ... calculate distance to seed
    dist = np.sqrt((X[pneighbors]-X[seed_ID])**2+(Y[pneighbors]-Y[seed_ID])**2)
    neighbors = pneighbors[dist<th_distance]

    # determine difference of orientation
    orient = np.zeros(len(neighbors))

    orient = abs(ON[neighbors] - seed_ON)
    orient[orient>90] = 180 - orient[orient>90]

    # create boolean based on thresholds --> label lineaments
    boolean = (abs(L[neighbors]-seed_LE)<th_length)& (orient < th_orient) & (neighbors!=seed_ID)
    labels = neighbors[boolean]
    dist = dist[dist<th_distance][boolean]
    RDs = np.array([])
    
    for i in labels:
        rd = calcRD(X[i], Y[i], seed_XY[0], seed_XY[1])
        RDs = np.append(RDs, rd)
    try:                    # get new seed
        new_seed = int(labels[dist==dist.max()])
        if new_seed == seed_ID: new_seed = None
    except:                 # no new seeds
        new_seed, new_seed2 = None, None
    
    # second seed
    try:
        dir_ns = RDs[dist == dist.max()][0]
        potential = labels[(RDs>=(dir_ns+90)) & (RDs<=(dir_ns+270))]
        potential = np.append(potential, labels[(RDs<=(dir_ns-90)) & (RDs>=(dir_ns-270))])
        pts = []
        for i in potential:
            pts.append(labels.tolist().index(i))
        new_seed2 = potential[dist[pts]==dist[pts].max()][0]
        if new_seed2 == seed_ID: new_seed2 = None
    except:
        new_seed2 = None
        
    return labels, new_seed, new_seed2    
    
def writeShapefile(lines, FIDs, IDs, label, seeds, sequence, source, previouslyClassified, L, OE, ON):
    shapeData = ogr.Open(lines, 1)
    layer = shapeData.GetLayer()
    FIDs_sorted = FIDs.copy()
    FIDs_sorted.sort()

    FIDs = IDs[FIDs_sorted.tolist()]
    i = 0   
        
    # required fields in shapefile
    newfields = ['label','order','seed','source', 'length', 'OE', 'ON'] #['id','label','order','seed','source']    
    for new_field in newfields:
        if layer.FindFieldIndex(new_field, 0) == -1:     # check if new fields already exist
            new_field = ogr.FieldDefn(new_field, ogr.OFTReal)   
            layer.CreateField(new_field)
    
    feature = layer.GetNextFeature()
    print '# lines classified with label %d: %d' % (label, len(FIDs))
    while feature:
        F = feature.GetFID() 
        if feature.GetField('label') == None:
            feature.SetField('label', 0)
                
        if previouslyClassified[F] == label:
            # remove earlier classification
            feature.SetField('source', 0)
            feature.SetField('label', 0)
            feature.SetField('order', 0)
            feature.SetField('seed', 0)
            feature.SetField('length', 0)        
            feature.SetField('OE', 0)        
            feature.SetField('ON', 0) 
        if i < len(FIDs):
            if F == FIDs[i]:
                feature.SetField('source', source[FIDs_sorted[i]])
                feature.SetField('label', label)
                feature.SetField('order', sequence[FIDs_sorted[i]])
                feature.SetField('seed', seeds[FIDs_sorted[i]])
                feature.SetField('length', L[FIDs_sorted[i]])        
                feature.SetField('OE', OE[FIDs_sorted[i]])        
                feature.SetField('ON', ON[FIDs_sorted[i]])        
                # Cleanup
                i +=1
        layer.SetFeature(feature)
        feature.Destroy()                   # clean cache
        feature = layer.GetNextFeature()
def checkShapefile(lines):
    shapeData = ogr.Open(lines, 1)
    layer = shapeData.GetLayer()
    feature = layer.GetNextFeature()
    while feature:
        print feature.GetFieldAsInteger('label')
        feature = layer.GetNextFeature()

def Clustre(IDs, x_start, x_end, y_start, y_end, X, Y, L, ON, OE, th_distance, th_orient, th_length, classified, plotYN='yes', seed_ID=0, label=1, output=None):
    
    # only use unclassified lines
    print 'seed ID: %d' % (seed_ID)
    IDs, x_start, x_end, y_start, y_end, X, Y, L, ON, OE = unclassified(IDs, x_start, x_end, y_start, y_end, X, Y, L, ON, OE, classified, label)
    
    try: seed_ID = np.argwhere(IDs == seed_ID)[0][0]
    except: print('\nWARNING: Seed already classified or not existing -> Unreliable results')
        
    FIDs = np.array([])
    seeds = [seed_ID]
    seeds_used = []
    
    FID = np.arange(len(IDs))
    order = 1
    sequence = np.zeros(len(IDs))
    source = np.zeros(len(IDs))

    while len(seeds)>0:
        # pick seed
        seed = seeds[0]             # pick seed
        seeds_used.append(seed)     # store seed in list of used seeds
        seeds = seeds[1:]           # remove seed from list of to-do seeds
        fid, new_seed, new_seed2 = labelling(seed, th_length, th_orient, th_distance, FIDs, FID, ON, OE, X, Y, L, order)
        sequence[fid] = order
        source[fid] = seed    
        order += 1
        for i in [new_seed, new_seed2]:
            if (i != None): seeds.append(i)
        for i in fid:
            if np.where(FIDs==i)[0].shape[0] < 1: FIDs = np.append(FIDs, i)


    if output != None:
        seeds = np.zeros(IDs.shape)
        seeds[seeds_used] = 1
        writeShapefile(output, FIDs, IDs, label, seeds, sequence, source, classified, L, OE, ON) 

    return FIDs.astype('int')
        
if __name__=='__main__':
    
    total   = len(sys.argv)

    if total != 8:      
      print "Invalid number of arguments"
      print "Define:\n lines=[location/to/shapefile.shp]\n th_distance=[distance threshold]\n th_orient=[threshold orientation]\n th_length=[threshold length]\n label=[classification]\n seed_ID=[initial seed FID]\n plotYN=['yes' or 'no' to (not) plot results]"
      sys.exit()
    
    for arg in sys.argv:
        if arg[-2:] != 'py':
            print arg
            exec(arg)

    # Read Shapefile
    box, IDs, x_start, x_end, y_start, y_end, X, Y, L, ON, OE, classified = readShapefile(lines)
    
    oldClassification = classified.copy()
    oldClassification[oldClassification == label] = 0   # old classification of label will be removed  
        
    # Run CLustre    
    classified = Clustre(IDs, x_start, x_end, y_start, y_end, X, Y, L, ON, OE, th_distance, th_orient, th_length, classified, seed_ID=seed_ID, label=label, output=lines, plotYN=plotYN)
    
    # PLot
    if plotYN=='yes':
        print 'Plotting...'
        plt.pause(0.5)
        
        fig = plt.figure(figsize=(10,10))
        for i in np.unique(oldClassification):
            i = i.astype('int')
            plt.plot((x_start[oldClassification==i],x_end[oldClassification==i]), (y_start[oldClassification==i],y_end[oldClassification==i]), c=np.random.rand(3,1), label=i)
        plt.plot((x_start[classified], x_end[classified]),(y_start[classified],y_end[classified]), color='r')
        plt.plot(X[seed_ID], Y[seed_ID],'*', markersize=10)
        
        plt.axis('equal')
        plt.xlim([x_start[x_start!=0].min(),x_start[x_start!=0].max()])
        plt.ylim([y_start[y_start!=0].min(),y_start[y_start!=0].max()])   
        title = 'seed ID: %d' % seed_ID
        plt.title(title)
        plt.savefig('preview.png',dpi=75)
    
    
    print 'Finished' 
    sys.exit(1)