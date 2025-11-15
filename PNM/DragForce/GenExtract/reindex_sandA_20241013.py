#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 22:59:37 2024

@author: htmt
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
import sys 
import os
from joblib import Parallel, delayed
sys.path.append('/home/e871/e871/mqu/GenExtract/')

from Pnm import pore_network_extraction as PNE
from Snm import solid_network_extraction as SNE
from ExtractionDNM import Mix_image,find_node_volume,find_interface,find_node_center

path_o='/home/e871/e871/mqu/Data/'
for k in ['SandA_N5000_3d','SandB_N5000_3d','SandC_N5000_3d']:
	path=path_o+k
	file_name='packing'
	
	data=np.fromfile(path+'/packing.raw',dtype=np.uint16)
	size=np.array(pd.read_csv(path+'/packing_size.dat',header=None)).flatten()
	print('Loading image data')
	resolution=size[0]
	size=size[-3:].astype(np.int16)
	data.shape=size
	data_list=np.unique(data)
	data_solid=np.copy(data)
	for i in tqdm(np.arange(len(data_list))):
		data_solid[data_solid==data_list[i]]=i
	data_list=np.vstack((data_list,np.arange(len(data_list)))).T
	data_list=pd.DataFrame(data_list)
	data_list.to_csv(path+'/index_list.csv')
	data_solid.astype(np.int32).tofile(path+'/packing_solid_solid.raw')
	data_pore=data.astype(bool)
	data.tofile(path+'/packing_o.raw')
	data_pore.tofile(path+'/packing.raw')
	
	
	PNE(path,file_name,resolution,size,path_abs='/home/e871/e871/mqu/GenExtract/',
		other_setting='0.05 0.98 0.7 0.5 1 1.1 4 0.15 1.75')
	path_pore  = path+'/'+file_name+'_pore_pore.raw'
	path_solid = path+'/'+file_name+'_solid_solid.raw'
	image,values,Solid_image,S_values,Pore_image,P_values=Mix_image(path_solid,path_pore,size)
	
	find_node_center(path,'solid',file_name,Solid_image,S_values)
	find_node_center(path,'pore',file_name,Pore_image,P_values)
	find_node_volume(path,'dual',file_name,image,values,core_number=96)
	find_interface(path,'dual',file_name,image,values,core_number=96)
	
	
	from scipy import ndimage
	import time
	import gc
	
	
	def Mix_image(path_solid,path_pore,size):
		picture_pore = np.fromfile(file=path_pore, dtype=np.int32)
		
		picture_solid = np.fromfile(file=path_solid, dtype=np.int32)
		
		picture_pore.shape = np.array(size)  #[1502,1502,1502]
		picture_solid.shape =np.array(size)#size
		
		array_pore = np.array(picture_pore)#[1:size[0]-1,1:size[1]-1,1:size[0]-1]
		array_solid = np.array(picture_solid)#[250:1250,250:1250,250:1250]
		del picture_pore,picture_solid
	
		#array_pore[array_pore<2]=0 #collate pore region, del useless region
	
		P_values = np.unique(array_pore)[1:] # obtain the No of pores 
	
		#values = values.tolist()
	
		max_pore = np.max(P_values)
	
		array_solid+=max_pore #recode the number of solid ball
		
		array_solid[array_solid<max_pore+1]=0 #collate region
		
		
		S_values=np.unique(array_solid)[1:] 
		
		values=np.concatenate((P_values,S_values)) # obtain the whole table
		#print(values)
		#array_pore=np.swapaxes(array_pore,2,0)
		array_mix = np.array(array_solid + array_pore, dtype=np.int32) #if max(values)>2**16 else np.array(array_solid + array_pore, dtype=np.int16)
	
		va=np.unique(array_mix)
		va=va[1:] if va[0]==0 else va
	
		if len(va)!=len(values):
			print('Error')
			sys.exit()
		else:
			print('No Error')
		#del array_pore,array_solid
		gc.collect()
		
		return array_mix,values,array_solid,S_values,array_pore,P_values
	t0 = time.time()
	
	
	path=path+'/'
	path_pore  = path+file_name+'_pore_pore.raw'
	path_solid = path+file_name+'_solid_solid.raw'
	
	index_data=pd.read_csv(path+'dual_network_interface_packing.csv')
	index_data=np.array(index_data)[:,:-1]
	solid_index_data=pd.read_csv(path+'solid_center_packing.csv')
	solid_index_data=np.array(solid_index_data)[:,2:-1]
	image,values,Solid_image,S_values,Pore_image,P_values=Mix_image(path_solid,path_pore,size)
	if os.path.exists(path+'throat_solid_pore')==False:    
		os.mkdir(path+'throat_solid_pore')
	
	def process_and_save(f, S_values, image, solid_index_data, index_data, path):  
		j = S_values[f]  
		image_c = np.copy(image)  
		image_c1 = np.copy(image)  
		image_c1[image_c1 != j] = 0  
		image_c1[image_c1 == j] = 1  
		structure1 = ndimage.generate_binary_structure(3, 1)  
		image_c1 = ndimage.binary_dilation(image_c1, structure=structure1)  
		image_c1 = image_c1.astype(int) * image_c  
		
		value = np.unique(image_c1)  
		
		for i in value[(value > 0) & (value < min(S_values))]:  
			structure1 = ndimage.generate_binary_structure(3, 1)  
			image_c2 = ndimage.binary_dilation(image_c1 == i, structure=structure1)  
			index = np.argwhere((image_c2 * image_c1).astype(bool))  
			data = np.ones_like(index) * solid_index_data[f]  
			index_df = pd.DataFrame(np.hstack((index, data)))  
			index_df.columns = ['x', 'y', 'z', 'solid_x', 'solid_y', 'solid_z']  
			index_throat = (index_data[(index_data[:, 1] == j) & (index_data[:, 2] == i)]).astype(int)  

			if len(index_throat) == 0:  
				continue  
			else:  
				index_throat = index_throat[0][0]  

			# Save CSV file with unique name  
			index_df.to_csv(f'{path}throat_solid_pore/{index_throat}__{j}.csv', index=False)  
	Parallel(n_jobs=128)(delayed(process_and_save)(f, S_values, image, solid_index_data, index_data, path) for f in tqdm(np.arange(len(S_values))))  

	'''
	for f in tqdm(np.arange(len(S_values))):
		j=S_values[f]
		image_c=np.copy(image)
		image_c1=np.copy(image)
		image_c1[image_c1!=j]=0
		image_c1[image_c1==j]=1
		structure1=ndimage.generate_binary_structure(3, 1)
		image_c1=ndimage.binary_dilation(image_c1,structure=structure1)
		image_c1=image_c1.astype(int)*image_c
		value=np.unique(image_c1)
		for i in  value[(value>0)&(value<min(S_values))]:
			structure1=ndimage.generate_binary_structure(3, 1)
			image_c2=ndimage.binary_dilation(image_c1==i,structure=structure1)
			#index=np.argwhere((image_c1==i)|(image_c2*image_c1).astype(bool))
			index=np.argwhere((image_c2*image_c1).astype(bool))
			data=np.ones_like(index)*solid_index_data[f]
			index=pd.DataFrame(np.hstack((index,data)))
			index.columns=['x','y','z','solid_x','solid_y','solid_z']
			index_throat=(index_data[(index_data[:,1]==j)&(index_data[:,2]==i)]).astype(int)
			if len(index_throat)==0:
				continue
			else:
				index_throat=index_throat[0][0]
			index.to_csv(path+'throat_solid_pore/'+'{}__{}.csv'.format(index_throat,j))
	'''
