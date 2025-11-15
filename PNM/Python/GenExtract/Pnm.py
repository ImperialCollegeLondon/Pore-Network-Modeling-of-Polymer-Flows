#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 09:43:07 2021

@author: htmt
"""

import os 
import fileinput as fi
import numpy as np
import shutil
import openpnm as op
import time

def pore_network_extraction(path,file_name,resolution,size,path_abs,direction='z',
							network_name=False,copy_phase=False,run=True,other_setting=[]):
	
	path_sample=path_abs+'/Sample.mhd'
	line=[]
	for l in fi.input(files=path_sample):
		line.append(l)
	fi.close()
	#line[5]='DimSize =    	1000	1000	1000\n'
	
	#path='/2/Documents/two_layers/generator_two_layers'
	
	data=np.fromfile(path+'/'+file_name+'.raw',dtype=np.uint8)
	data.shape=size
	#data=data[100:1100,100:1100,100:1100]
	#resolution=50e-6
	network_name=file_name+'_pore'
	size=data.shape
	data[data>0]=1
	if copy_phase:
		
		data.tofile(path+'/'+network_name+'.raw')
	if 'z' in direction:
		picture=np.swapaxes(data,2,0)    
	elif 'y' in direction:
		picture=np.swapaxes(data,2,1) 
	elif 'x' in direction:
		picture=picture
	picture.tofile(path+'/'+network_name+'_inv.raw')
	del picture, data
	path+='/pore_network'
	if os.path.exists(path)==False:    
		os.mkdir(path)
	
	line[5]='DimSize =  	'+str(size[2])+'	'+str(size[1])+'	'+str(size[0])+'\n'
	line[6]='ElementSize = 	'+str(resolution)+' 	'+str(resolution)+' 	'+str(resolution)+' \n'
	
	line[9]='ElementDataFile = ../'+file_name+'.raw\n'
	line[11]='DefaultImageFormat = .raw\n'
	line[13]='direction '+direction+'\n'
	line[15]='pore 0 0\n'
	if len(other_setting)!=0:
		line.append('medialSurfaceSettings:'+other_setting+'\n')
		
		
	if os.path.exists(path+'/Sample.mhd')==False:    
		sample_w=open(path+'/Sample.mhd','a+')
	else:
		os.remove(path+'/Sample.mhd')
		sample_w=open(path+'/Sample.mhd','a+')
	for j in line:
		sample_w.write(j)
	sample_w.close()
	a = os.getcwd()    
	
	## Address of the mhd file for image
	mhd_Address = path+'/Sample.mhd'
	## Address of the network extraction code 
	NE_Address = path_abs+"/pnextract"
	if os.path.exists(mhd_Address+'/pnextract')==False:
		shutil.copyfile(NE_Address,os.path.dirname(mhd_Address)+'/pnextract')
	else:
		print(f'{mhd_Address}/pnextract exists')
	os.chmod(os.path.dirname(mhd_Address)+'/pnextract', 0o775)
	os.chdir(os.path.dirname(mhd_Address))
	if run:
		cmd = os.system("./pnextract "+ os.path.basename(mhd_Address))    
	if os.path.exists(a+os.path.dirname(mhd_Address)[1:]+'/pnextract')==True:
		os.remove(a+os.path.dirname(mhd_Address)[1:]+'/pnextract') 
	elif os.path.exists(os.path.dirname(mhd_Address)+'/pnextract')==True:
		os.remove(os.path.dirname(mhd_Address)+'/pnextract') 
	os.chdir(a)
	
	path = os.path.dirname(mhd_Address)#'./pore_network'#Path('../_fixtures/ICL-Sandstone(Berea)/')
	#project = op.io.Statoil.load(path=os.path.dirname(mhd_Address), prefix=file_name)
	project = op.io.network_from_statoil(path=os.path.dirname(mhd_Address), prefix=file_name)
	pn = project.network
	pn.name = 'pore'
	#size+=2 # please note the order has change after pnextract[Num,High,Width]
	
	path_pore=path+'/'+file_name+'_VElems.raw'
	picture_pore = np.fromfile(file=path_pore, dtype=np.int32)
	size=np.array(size)
	if 'z' in direction:
		size=size[[2,1,0]]    
	elif 'y' in direction:
		size=size[[0,2,1]] 
	elif 'x' in direction:
		size=size
	picture_pore.shape = np.array(size)+2      #[1502,1502,1502]
	array_pore = np.array(picture_pore)[1:size[0]+1,1:size[1]+1,1:size[2]+1]-1
	array_pore[array_pore<1]=0
	pore_number=np.unique(array_pore)
	if 'z' in direction:
		
		array_pore=np.swapaxes(array_pore,2,0)   
	elif 'y' in direction:
		array_pore=np.swapaxes(array_pore,2,1)
	elif 'x' in direction:
		array_pore=array_pore
	
	pore_index=list(set(np.arange(max(pore_number))+1)-set(pore_number))
	if len(pore_index)>0:
		print('it needs to be modified')
		for i in pore_index:
			index=pn['pore.coords'][i]/float(resolution)
			array_pore[int(index[0]),int(index[1]),int(index[2])]=i
	else:
		print('it does not need to be modified')    
	 
	array_pore.tofile(path+'/../'+network_name+'_pore.raw')
	
	print('finish the part')
