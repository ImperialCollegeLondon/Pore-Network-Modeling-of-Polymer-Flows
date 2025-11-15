#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 01:11:41 2024

@author: htmt
"""
import sys 
sys.path.append('/home/e871/e871/mqu/GenExtract/')
import argparse 

from Pnm import pore_network_extraction as PNE
from Snm import solid_network_extraction as SNE
from ExtractionDNM import Mix_image,find_node_volume,find_interface,find_node_center

parser=argparse.ArgumentParser(description='Dual netowrk extraction')
parser.add_argument('-p','--path',default='./')
parser.add_argument('-n','--name',default='sphere')
parser.add_argument('-s','--size',type=int,default=[1000,500,500], help='input batch size')
parser.add_argument('-r','--resolution',default='1e-6')
args=parser.parse_args()
print(args)
path=args.path
file_name=args.name
size=args.size
resolution=args.resolution
PNE(path,file_name,resolution,size,path_abs='/home/e871/e871/mqu/GenExtract/')
SNE(path,file_name,resolution,size)

path_pore  = path+file_name+'_pore_pore.raw'
path_solid = path+file_name+'_solid_solid.raw'
image,values,Solid_image,S_values,Pore_image,P_values=Mix_image(path_solid,path_pore,size)

find_node_center('solid',file_name,Solid_image,S_values)
find_node_center('pore',file_name,Pore_image,P_values)
find_node_volume('dual',file_name,image,values)
find_interface('dual',file_name,image,values)