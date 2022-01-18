%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A wavenumber-domain iterative approach for 3D imaging of magnetic anomalies and gradients
% Author: Yatong Cui, Lianghui Guo (guolh@cugb.edu.cn)
% Organization: China University of Geosciences (Beijing), School of Geophysics and Information Technology
% Compiled version: MATLAB R2017b
% Reference:
%       Cui Y T, Guo L H. A wavenumber-domain iterative approach for 3D imaging of magnetic
%       anomalies and gradients with depth constraints. Journal of Geophysics and Engineering, 
%       2019, 16(6): 1032-1047.
% Description of the input parameters: 
%       infile_ano: observed anomaly, unit: nT
%       inclination, declination£ºgeomagnetic inclination and declination 
%                  The direction of magnetization is consistent with the direction of geomagnetic field
%       n: power of vertical derivative, n¡Ê[1,10]
%       N: depth scaling factor£¬N¡Ý1
%       zmin, zmax: range of imaging depth direction
%       dz: spacing of depth direction
%       iter: number of iterations
%       m1, m2: physical property range constraint, unit: A/m
% Description of the output parameters: 
%       outfile_msh: imaging model mesh file
%       outfile_mod: imaging model property file, unit: A/m
%       outfile_rms: root mean square file, unit: nT
%       outfile_res: residual anomaly file, unit: nT
%       outfile_for: calculated anomaly file, unit: nT
% Description of primary identifiers£º
%       x, y: x, y verctor
%       nx, ny: number of points in x and y directions
%       dx, dy: spacing in x and y directions
%       npts: extension points
%       u0: vacuum magnetic permeability, unit: Henry/m
%       F£ºgeomagnetic field direction vector
%       M£ºmagnetization direction vector
%       kx, ky, k: wavenumber
%       obs: observed anomaly
%       cal: calculated anomaly
%       mcal: imaging model
%       dm: magnetization disturbance
%       res: residual anomaly
%       error: root mean square
% Description of subroutine function: 
%       readgrd.m: read surfer text grd file
%       wave2d.m: calculate wavenumber
%       imaging_Ut.m: imaging
%       forward_Ut.m: forward
%       extend_copy2d.m: copy edge extension
%       property.m: physical property range constraint
%       savebln.m: save rms file
%       savegrd.m: save surfer text grd file
%       savemod.m: save model file
%       savemsh.m: save mesh file      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%%%%%%%%%%% I/O parameters %%%%%%%%%%%%%
infile_ano = 'Ut_2D.grd'; 
inclination = 60;
declination = 45;
n = 4;
N = 2.5; 
zmin = 0; zmax = 2000; 
dz = 50; 
iter = 5;
m1 = 0; m2 = 1; 
outfile_msh = 'img_model.msh';
outfile_mod = 'img_model.mod'; 
outfile_rms = 'img_model_rms.bln'; 
outfile_res = 'img_model_res.grd'; 
outfile_for = 'img_model_for.grd'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[obs,x,y,nx,ny,dx,dy] = readgrd(infile_ano);
nz = floor((zmax-zmin)/dz+1); z = zmin:dz:zmax;
nmax = max([nx ny]);
npts = 2^nextpow2(nmax);
[F,M] = getINDE(inclination,declination);
u0 = 4*pi*10^(-7);error = zeros(iter,1);
% imaging
[k,kx,ky] = wave2d(npts,dx,dy);
mcal = imaging_Ut(obs,u0,M,F,n,N,nx,ny,nz,dz,z,npts,k,kx,ky); 
mcal = property(mcal,nx,ny,nz,m1,m2); 
cal = forward_Ut(mcal,u0,M,F,nx,ny,nz,dz,z,npts,k,kx,ky);
res = obs-cal;
for i = 1:iter
    dm = imaging_Ut(res,u0,M,F,n,N,nx,ny,nz,dz,z,npts,k,kx,ky); 
    mcal = mcal+dm;
    mcal = property(mcal,nx,ny,nz,m1,m2);
    cal = forward_Ut(mcal,u0,M,F,nx,ny,nz,dz,z,npts,k,kx,ky);
    res = obs-cal; 
    error(i) = rms(res(:));
 end
% save
savemsh(x,y,z,dx,dy,dz,nx,ny,nz,outfile_msh)
savemod(mcal,nx,ny,nz,outfile_mod);
savebln(error,outfile_rms);  
savegrd(res,x,y,nx,ny,outfile_res);
savegrd(cal,x,y,nx,ny,outfile_for);