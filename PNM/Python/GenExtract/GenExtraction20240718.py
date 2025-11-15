#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 01:11:41 2024

@author: htmt
"""
import sys 
sys.path.append('/home/htmt/Documents/GenExtract/')
import argparse 
import numpy as np
from Pnm import pore_network_extraction as PNE
from Snm import solid_network_extraction as SNE
from ExtractionDNM import Mix_image,find_node_volume,find_interface,find_node_center
from scipy import ndimage
from joblib import Parallel, delayed
from tqdm import tqdm
import gc
from skimage import measure
import os 
import fileinput as fi
import shutil
import openpnm as op
import numpy as np
from scipy import ndimage as ndi
from skimage import feature
from skimage.segmentation import watershed
import time

parser=argparse.ArgumentParser(description='Dual netowrk extraction')
parser.add_argument('-p','--path',default='./OrderedPattern_ppp50')
parser.add_argument('-n','--name',default='packing')
parser.add_argument('-s','--size',metavar='N',type=int,default=[330,330,330],nargs='+', help='input batch size')
parser.add_argument('-r','--resolution',default='1e-5')
parser.add_argument('-d','--diameter',type=int,default=25)
parser.add_argument('-c','--condition',default='pnextract')
args=parser.parse_args()
print(args)
path=args.path
file_name=args.name
size=np.array(args.size) #eval(''.join(args.size)))
resolution=args.resolution
if 'pnextract' in args.condition:
    PNE(path,file_name,resolution,size,path_abs='/home/htmt/Documents/GenExtract/')
elif 'genextract' in args.condition:
    SNE(path,file_name,resolution,size,network_type='pore',direction='x',distance_seed=args.diameter)
'''
SNE(path,file_name,resolution,size,distance_seed=args.diameter)
'''
path_pore  = path+'/'+file_name+'_pore_pore.raw'
path_solid = path+'/'+file_name+'_solid_solid.raw'
image,values,Solid_image,S_values,Pore_image,P_values=Mix_image(path_solid,path_pore,size)

find_node_center(path,'solid',file_name,Solid_image,S_values)
find_node_center(path,'pore',file_name,Pore_image,P_values)
find_node_volume(path,'dual',file_name,image,values)
find_interface(path,'dual',file_name,image,values)

