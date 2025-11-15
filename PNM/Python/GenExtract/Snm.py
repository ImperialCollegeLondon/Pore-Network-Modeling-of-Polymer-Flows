#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 15:09:04 2021

@author: mingliangqu
"""

import numpy as np
from scipy import ndimage as ndi
from skimage import feature
from skimage.segmentation import watershed
import time

def solid_network_extraction(path,file_name,resolution,size,network_type='solid',network_name=False,distance_seed=5,copy_phase=False):
    t0 = time.time()
    #path='./'
    #file_name='sphere_stacking_500_500_2000'
    network_name= network_name if network_name else file_name+'_'+network_type
    #resolution=7.87e-6
    #size=np.array([2000,500,500])
    
    image=np.fromfile(path+'/'+file_name+'.raw',dtype=np.uint8)
    image.shape=size
    #image=image[100:1100,100:1100,100:1100]
    if copy_phase:
        image.tofile(path+'/'+network_name+'.raw')
    
    
    #picture=np.swapaxes(image,2,0)    
    #picture.tofile(path+'/'+network_name+'_inv.raw')
    #del picture
    
    tend = time.time()
    print("time cost：%.6fs" % (tend - t0))
    t0 = time.time()
    #image[image<-1]=0
    #image=np.abs(image)
    if network_type == 'pore':
        image=(~image.astype(bool)).astype(np.uint8)
    distance = ndi.distance_transform_edt(image) #distance_map
    local_maxi=feature.peak_local_max(distance,min_distance=distance_seed, 
                                      indices=False, 
                                      exclude_border=False,
                                      footprint=np.ones((3, 3, 3)),
                                      labels=image)   #seek peak
    #plt.imshow(ndi.binary_dilation(local_maxi,np.ones((3,3,3)),iterations=5)[:,2,:]) visualization
    markers = ndi.label(local_maxi)[0] #markers
    labels =watershed(image=-distance, markers=markers, #compactness=0.1,
                      mask=image) #Watershed algarithm
    
    labels.tofile(path+'/'+network_name+'_'+network_type+'.raw')
    
    tend = time.time()
    print("time cost：%.6fs" % (tend - t0))
