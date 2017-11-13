function [output] = load2decomp_RAV(output,folder,fs3D,fs2D,nx,ny,nz,outsteps)
% =========================================================================
% load running average (RAV) files
% 
% CALL:  load2decomp_RAV(output,folder,fs3D,fs2D,nx,ny,nz,outsteps)
% INPUT: output - INOUT data structre
%        fs3D - variable in the 3D RAV files
%        fs2D - variabel in the 2D RAV files
%        folder - path to data
%        nx,ny,nz - grid resolution
%        outsteps - array of output steps for RAV files
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================

len3D=nx*ny*nz;
len2D=nx*ny;
n=length(outsteps);
nt_total=0;

for kk=1:numel(fs3D)
    output.data3D.(fs3D{kk})=[];
end
for kk=1:numel(fs2D)
    output.data2D.(fs2D{kk})=[];
end

for ii=1:n
    % read data info file
    fn=sprintf('%s/RAV_%09i_info.txt',folder,outsteps(ii));
    infodata=load(fn);
    if(infodata(1)~=nx ||infodata(2)~=ny || infodata(3)~=nz)
        return
    end
    nt = infodata(4);
    nt_total=nt_total+nt;
    
    % read 3D data file
    fn=sprintf('%s/RAV_%09i_3D.dat',folder,outsteps(ii));
    fid = fopen(fn, 'rt');
    tmp = [fread(fid,'double')];
    fclose(fid);
    
    displ=0;
    % reshape data and save to structure
    for kk=1:numel(fs3D)
        output.data3D.(fs3D{kk})(:,:,:,ii)=reshape(tmp(displ+1:displ+len3D),nx,ny,nz);
        displ=displ+len3D;
    end
    
    % read 2D data file
    fn=sprintf('%s/RAV_%09i_2D.dat',folder,outsteps(ii));
    fid = fopen(fn, 'rt');
    tmp = [fread(fid,'double')];
    fclose(fid);
    
    displ=0;
    % reshape data and save to structure
    for kk=1:numel(fs2D)
        output.data2D.(fs2D{kk})(:,:,ii)=reshape(tmp(displ+1:displ+len2D),nx,ny);
        displ=displ+len2D;
    end
end
output.nt=nt_total;
end
