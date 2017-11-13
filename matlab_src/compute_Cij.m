function sout = compute_Cij(s)
% =========================================================================
% compute dispersive stresses 
%
% INPUT: s - data structure containing u,v,w 
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
C11=zeros(size(s.u,3),size(s.u,4));
C12=zeros(size(s.u,3),size(s.u,4));
C22=zeros(size(s.u,3),size(s.u,4));

C13=zeros(size(s.w,3),size(s.w,4));
C23=zeros(size(s.w,3),size(s.w,4));
C33=zeros(size(s.w,3),size(s.w,4));

u = colocate_var(s.u,'w');
v = colocate_var(s.v,'w');

for t=1:size(s.u,4)
    C11(:,t)=get_stress(s.u(:,:,:,t),s.u(:,:,:,t));
    C12(:,t)=get_stress(s.u(:,:,:,t),s.v(:,:,:,t));
    C22(:,t)=get_stress(s.v(:,:,:,t),s.v(:,:,:,t));
        
    C13(:,t)=get_stress(u(:,:,:,t),s.w(:,:,:,t));
    C23(:,t)=get_stress(v(:,:,:,t),s.w(:,:,:,t));
    C33(:,t)=get_stress(s.w(:,:,:,t),s.w(:,:,:,t));
end

sout=struct('C11',C11,'C12',C12,'C22',C22,'C13',C13,'C23',C23,'C33',C33);

end

function out = get_stress(ui,uj)
out=squeeze(mean(mean(ui.*uj,1),2)) -...
    squeeze(mean(mean(ui,1),2)).*squeeze(mean(mean(uj,1),2));
end