function [Fp_x, Fp_y, Fp_z, totalArea] = integrateSurfacePressure(porePressure, xp, yp, zp, rp, x_voxel, y_voxel, z_voxel, dx, colorvec)
dist = vecnorm([x_voxel-xp, y_voxel-yp, z_voxel-zp], 2, 2);
Fp_x = sum((xp-x_voxel)./dist)*porePressure*dx*dx;
Fp_y = sum((yp-y_voxel)./dist)*porePressure*dx*dx;
Fp_z = sum((zp-z_voxel)./dist)*porePressure*dx*dx;
totalArea = length(x_voxel)*dx*dx;






