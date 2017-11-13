% =========================================================================
% Example of matlab script, loading running average data from 2decomp LES 
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================

% =========================================================================
set(0,'DefaultAxesFontSize',10)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
% =========================================================================
ustar=0.45;z_i=1000;
z0=0.1/z_i;
dt=0.2;
nx=128;ny=128;nz=128;
lx=2*pi;ly=lx;lz=1;
dx=lx/nx;dy=ly/ny;dz=lz/nz;
x=0:dx:lx-dx;y=0:dy:ly-dy;
z_w=0:dz:lz-dz;z_uvp=z_w+0.5*dz;
z_th=0:0.001:lz;u_prof_th=1/0.4*log(z_th/z0);

% outsteps - array of output steps for RAV files
stat_start=1.5E5;stat_end=2E5;stat_np=2.5E4;
outsteps = stat_start+stat_np:stat_np:stat_end;


data=struct();
path='/scratch/general/lustre/u0917142/debug/';
folder = [path,'/debug_test2'];cf='test';
data.(cf)=simulation_stats_check(nx,ny,nz,lx,ly,lz,dt,z0,ustar,z_i,stat_start,stat_end,stat_np,folder,cf,true);