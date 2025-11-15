#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 20:52:49 2021

@author: htmt
"""
from scipy import ndimage
import numpy as np
import pandas as pd
from tqdm import tqdm
import time
import sys
#from skimage import feature
from joblib import Parallel, delayed
import gc
from skimage import measure
'''
def gradient(image):
    grd_z, grd_y, grd_x = np.gradient(image,edge_order=1)
    grd_z, grd_y, grd_x = grd_z.astype(np.float32), grd_y.astype(np.float32), grd_x.astype(np.float32) 
    grad = np.abs(grd_z) + np.abs(grd_y) + np.abs(grd_x)
    return grad
'''
def find(image, i):
    region_i=np.argwhere(image == i)
    shape=image.shape
    if region_i.size != 0 :

        zmax = np.max(region_i[:, 0]) + 2 #if np.max(region_i[:, 0]) +1<shape[0] else shape[0]
        zmin = np.min(region_i[:, 0]) - 1 #if np.min(region_i[:, 0]) -1>0 else 0
        xmax = np.max(region_i[:, 1]) + 2 #if np.max(region_i[:, 1]) +1<shape[1] else shape[1]
        xmin = np.min(region_i[:, 1]) - 1 #if np.min(region_i[:, 1]) -1>0 else 0
        ymax = np.max(region_i[:, 2]) + 2 #if np.max(region_i[:, 2]) +1<shape[2] else shape[2]
        ymin = np.min(region_i[:, 2]) - 1 #if np.min(region_i[:, 2]) -1>0 else 0
        index=np.array([zmin,xmin,ymin])
        index[index<0]=0
        index1=np.array([zmax,xmax,ymax])
        index1=(index1>shape)*(np.array(shape)+1)+(index1<=shape)*index1
        region=np.copy(image[index[0]:index1[0],index[1]:index1[1],index[2]:index1[2]])
        #region = np.pad(region, ((1, 1), (1, 1), (1, 1)), 'constant', constant_values=(-1, -1)) #pad region should be discussed
        '''
        if ~np.any(~(np.array(np.shape(region))-1>0)):
            grad=gradient(region)
            tempgz,tempgx,tempgy  = np.gradient(region,edge_order=1)
            #gz,gy,gx =diff_interface(region)
                
        else:
            return 0
        '''
        return region, index#,grad, tempgz, tempgy, tempgx, #,gz,gy,gx
    else:
        return 0
def calc_surf(image,i):
    outcome0 = find(image, i)
    if outcome0 != 0:
        outcome1 = calculate_interface_surface(outcome0,i)
        return outcome1
    else:
        return np.array([i, 0, 0])
    
def calculate_interface_surface(outcome0,i):
    region=np.copy(outcome0[0])
    verts,faces,normals,values=measure.marching_cubes((region==i))
    area=measure.mesh_surface_area(verts,faces)
    gc.collect()
    region_0=np.copy(region)
    region_0[region_0!=i]=0
    region[region==i]=0
    structure1=ndimage.generate_binary_structure(3, 3)
    region_0=ndimage.binary_dilation(region_0,structure=structure1)
    region=region_0*region
    values=np.unique(region)
    if len(values[values>0])<=1:
        table=[i,-1,-1]
        
    else:
        table = Parallel(n_jobs=2)(delayed(measure_surf_test)(np.copy(outcome0[0]),i,j) for j in values[values>0])
        table=np.array(table)
        table[:,1]=(table[:,1]+area-table[:,2])/2
        table=table[:,:2]
        
        del outcome0
        del region
        gc.collect()
        table=np.c_[np.ones(len(values[values>0]))*i,table].astype(np.float32)
        table=table[table[:,2]>0]
        #table[table[:,2]<0]=
        if len(table)<=1:
            table=[i,-1,-1]
    
    return table
def measure_surf_test(region,i,j):
    verts1,faces1,normals1,values1=measure.marching_cubes((region==j))
    area1=measure.mesh_surface_area(verts1,faces1)
    verts2,faces2,normals2,values2=measure.marching_cubes((region==i)+(region==j))
    area2=measure.mesh_surface_area(verts2,faces2)
    
    return np.array([j,area1,area2])
    
def measure_surf(image,i):
    outcome = find(image, i)
    if outcome != 0:
        outcome1 = calculate_surface(outcome,i)
        outcome1[1]/=2
        return outcome1
    else:
        return np.array([i, 0])

def calculate_surface(array,i):
    region=array[0]
    '''
    distance=((region_index-center_0[0]).dot(norm_cen)).astype(np.int32)
    shadow_index=((region_index-center_0[0])-distance.reshape([len(distance),1])*norm_cen/np.linalg.norm(norm_cen)**2).astype(np.int32)
    '''
    
    del array
    verts,faces,normals,values=measure.marching_cubes(region==i)
    area=measure.mesh_surface_area(verts,faces)
    
    return [i,area]
  

def cal_volume(array,i):
    region=np.copy(array[0])
    del array
    Bin_region=region/i
    Bin_region[Bin_region!=1]=0

    volume=np.count_nonzero(Bin_region)
    

    verts,faces,normals,values=measure.marching_cubes(Bin_region==1)
    area=measure.mesh_surface_area(verts,faces)
    #print(total)
    return volume,area 

def calc_volum(image, i):
    outcome = find(image, i)
    if outcome != 0:
        volume_surf = cal_volume(outcome,i)
        #print(i) 
        return [i,volume_surf[0],volume_surf[1]]
    else:
        return [i,0,0]

def test(region,i):
    import pyvista as pv
    verts,faces,normals,values=measure.marching_cubes(region==i)
    faces=np.column_stack((np.ones(len(faces),)*3,faces)).astype(int)
    grid=pv.PolyData(verts,faces)
    grid['Normals']=normals
    grid['values']=values
    plotter=pv.Plotter()
    plotter.add_mesh(grid,scalars='values')
    plotter.show()

def center(image, i):
    #region_0=np.copy(region)
    region=np.copy(image[0])
    index=np.copy(image[1])
    del image
    region[region!=i]=0
    region=region/i
    distance    = ndimage.distance_transform_edt(region).astype(np.float16)
    del region
    centerindex = np.average(np.argwhere(distance==np.max(distance)),axis=0)
    radius = np.array([distance[int(centerindex[0]),int(centerindex[1]),int(centerindex[2])]])
    if radius==0:
        radius=1
    
    centerindex += index
    result = np.append(centerindex, radius)
    result = np.append(i,result)
    gc.collect()
    return result

def find_center(image,i):
    outcome0 = find(image, i)
    if outcome0 != 0:
        outcome1 = center(outcome0,i)
        return outcome1
    else:
        return np.array([i, 0, 0, 0, 0])


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




def find_node_center(path,net_type,name,image,index=[],core_number=24):
    t0 = time.time()
    gc.collect()
    index=np.unique(image) if len(index)==0 else index
    arrays = Parallel(n_jobs=core_number,prefer='threads')(delayed(find_center)(image,i) for i in tqdm(index))
    tend = time.time()
    #cen_table=np.vstack(arrays)
    cen_table=pd.DataFrame(np.vstack(arrays))
    
    cen_table.to_csv(path+'/'+net_type+'_center_'+name+'.csv')
    print('running time of finding '+name+' center：%.6fs' % (tend - t0))
def find_node_volume(path,net_type,name,image,index=[],core_number=24):
    t0 = time.time()    
    gc.collect()   
    index=np.unique(image) if len(index)==0 else index
    result_v= Parallel(n_jobs=core_number,prefer='threads')(delayed(calc_volum)(image, i) for i in tqdm(index[index>0]))
    result_v=np.vstack(result_v).astype(np.int32)
    
    result_v=pd.DataFrame(result_v)
    
    result_v.to_csv(path+'/'+net_type+'_network_volume_'+name+'.csv')
    tend = time.time()
    print('running time of calculation volume：%.6fs' % (tend - t0))
    
def find_interface(path,net_type,name,image,index=[],core_number=24):
    t0 = time.time()    
    gc.collect() 
    index=np.unique(image) if len(index)==0 else index
    arrays = Parallel(n_jobs=core_number,prefer='threads')(delayed(calc_surf)(image ,i) for i in tqdm(index[index>0]))
    
        
    table=np.vstack(arrays).astype(np.float32)
    table=table[(table[:,2]>0)&(table[:,1]>0)]
    
    table1=np.vstack((table,table[:,[1,0,2]]))
    table1=table1[table1[:,1].argsort()]
    table1=table1[np.lexsort((table1[:,0],))]  
    table2=pd.DataFrame(table1[:,:2])
    table2=table2.apply(tuple,axis=1)
    table2=pd.DataFrame(table2.values.T)
    table2[1]=table1[:,2]
    table2=table2.drop_duplicates(subset=0)
    table=pd.DataFrame(table2[0].tolist(),index=table2.index)
    table[2]=table2[1]
    table=table.reset_index()
    table=table.drop('index',axis=1,inplace=False)
    table=pd.DataFrame(table)
    table[0]=table[0].astype(np.int32)
    table[1]=table[1].astype(np.int32)
    table.columns=['index_1','index_2','Area']
    table.to_csv(path+'/'+net_type+'_network_interface_'+name+'.csv')
    
    
    
    tend = time.time()
    
    print("running time of calculation interface：%.6fs" % (tend - t0))    
    
'''
path='./'
file_name='sphere_stacking_500_500_2000'
path_pore  = path+file_name+'_pore_pore.raw'
path_solid = path+file_name+'_solid_solid.raw'
size=[2000,500,500]
image,values,Solid_image,S_values,Pore_image,P_values=Mix_image(path_solid,path_pore,size)
#Pore_image.tofile('./bead_stacking_500_12_uint8_pore_c.raw')

#P_values = np.unique(image)
name=file_name
find_node_center('solid',name,Solid_image,S_values)
find_node_center('pore',name,Pore_image,P_values)
find_node_volume('dual',name,image,values)
find_interface('dual',name,image,values)
'''

