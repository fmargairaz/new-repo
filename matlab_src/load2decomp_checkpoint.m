function chpt = load2decomp_checkpoint(folder,nx,ny,nz,outstep,LSD_flag)
% =========================================================================
% load checkpoint files
% 
% CALL:  load2decomp_checkpoint(folder,nx,ny,nz,outstep,LSD_flag)
% INPUT: folder - path to data
%        nx,ny,nz - grid resolution
%        outstep - output step for checkpoint files
%        LSD_flag - flag for Lagrangian scale dependent variable
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================

fs3D1={'u','v','w','RHSx','RHSy','RHSz'};

if(LSD_flag)
    fs3D2={'Cs_opt2','F_LM','F_MM','F_QN','F_NN'};
else
    fs3D2={'Cs_opt2'};
end

% read 3D data file
fn=sprintf('%s/checkpoint_%09i',folder,outstep);
fprintf('==============================================================\n')
fprintf('reading checkpoint file:\n%s\n',fn)
fprintf('==============================================================\n')

fid = fopen(fn, 'rt');
tmp = [fread(fid,'double')];
fclose(fid);

len3D=nx*ny*(nz+1);
displ=0;
% reshape data and save to structure
for kk=1:numel(fs3D1)
    chpt.(fs3D1{kk})(:,:,:)=reshape(tmp(displ+1:displ+len3D),nx,ny,nz+1);
    displ=displ+len3D;
end

len3D=nx*ny*(nz+1);
% reshape data and save to structure
for kk=1:numel(fs3D2)
    chpt.(fs3D2{kk})(:,:,:)=reshape(tmp(displ+1:displ+len3D),nx,ny,nz+1);
    displ=displ+len3D;
end

end