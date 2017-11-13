function [output] = get_data_struct_RAV(folder,nx,ny,nz,outsteps)
% =========================================================================
% load data3D and data2D structure from running averages
% 
% CALL:  get_data_struct_RAV(folder,nx,ny,nz,outsteps)
% INPUT: folder - path to data
%        nx,ny,nz - grid resolution
%        outsteps - array of output steps for RAV files
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
output = struct('data3D',[],'data2D',[]);

fs3D1 = {'u','v','w','p',...
    'u2','v2','w2','u3','v3','w3','u4','v4','w4',...
    'uv','uw','vw','txx','tyy','tzz','txy','txz','tyz',...
    'dudz','dvdz','nut','cs'};
fs2D1={'ustar'};
output=load2decomp_RAV(output,folder,fs3D1,fs2D1,nx,ny,nz,outsteps);

fs={fs3D1{:}};
for kk=1:numel(fs)
    output.data3D.(fs{kk})=squeeze(sum(output.data3D.(fs{kk}),4))/output.nt;
end

fs={fs2D1{:}};
for kk=1:numel(fs)
    output.data2D.(fs{kk})=squeeze(sum(output.data2D.(fs{kk}),3))/output.nt;
end

fs_tij = {'txx','tyy','tzz','txy','txz','tyz'};
for kk=1:numel(fs_tij)
    output.data3D.(fs_tij{kk})=-output.data3D.(fs_tij{kk});
end

output.data3D.w_uvp=0.5*(output.data3D.w(:,:,1:nz-1)+output.data3D.w(:,:,2:nz));
output.data3D.w_uvp(:,:,nz) = 0.5*output.data3D.w(:,:,nz-1);
output.data3D.w2_uvp=0.5*(output.data3D.w2(:,:,1:nz-1)+output.data3D.w2(:,:,2:nz));
output.data3D.w2_uvp(:,:,nz) = 0.5*output.data3D.w2(:,:,nz-1);
output.data3D.w3_uvp=0.5*(output.data3D.w3(:,:,1:nz-1)+output.data3D.w3(:,:,2:nz));
output.data3D.w3_uvp(:,:,nz) = 0.5*output.data3D.w3(:,:,nz-1);
output.data3D.w4_uvp=0.5*(output.data3D.w4(:,:,1:nz-1)+output.data3D.w4(:,:,2:nz));
output.data3D.w4_uvp(:,:,nz) = 0.5*output.data3D.w4(:,:,nz-1);

end