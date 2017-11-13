function initcond = load2decomp_initialcond(folder,nx,ny,nz)
% =========================================================================
% load initial condition files
% 
% CALL:  load2decomp_initialcond(folder,nx,ny,nz)
% INPUT: folder - path to data
%        nx,ny,nz - grid resolution
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================

fs3D1={'u','v','w'};

% read 3D data file
fn=sprintf('%s/initial_conditions.dat',folder);
fprintf('==============================================================\n')
fprintf('reading initial conditions file:\n%s\n',fn)
fprintf('==============================================================\n')

fid = fopen(fn, 'rt');
tmp = [fread(fid,'double')];
fclose(fid);

len3D=nx*ny*nz;
displ=0;
% reshape data and save to structure
for kk=1:numel(fs3D1)
    initcond.(fs3D1{kk})(:,:,:)=reshape(tmp(displ+1:displ+len3D),nx,ny,nz);
    displ=displ+len3D;
end

end