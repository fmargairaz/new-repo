function sout = compute_Rij(s)
% =========================================================================
% compute resolved stresses 
%
% INPUT: s - data structure containing u,v,w 
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================

R11=zeros(size(s.u));
R12=zeros(size(s.u));
R22=zeros(size(s.u));

R13=zeros(size(s.w));
R23=zeros(size(s.w));
R33=zeros(size(s.w));

u = colocate_var(s.u,'w');
v = colocate_var(s.v,'w');

for t=1:size(s.u,4)
    R11(:,:,:,t)=s.u2(:,:,:,t)-s.u(:,:,:,t).*s.u(:,:,:,t);
    R12(:,:,:,t)=s.uv(:,:,:,t)-s.u(:,:,:,t).*s.v(:,:,:,t);
    R22(:,:,:,t)=s.v2(:,:,:,t)-s.v(:,:,:,t).*s.v(:,:,:,t);
        
    R13(:,:,:,t)=s.uw(:,:,:,t)-u(:,:,:,t).*s.w(:,:,:,t);
    R23(:,:,:,t)=s.vw(:,:,:,t)-v(:,:,:,t).*s.w(:,:,:,t);
    R33(:,:,:,t)=s.w2(:,:,:,t)-s.w(:,:,:,t).*s.w(:,:,:,t);
end

sout=struct('R11',R11,'R12',R12,'R22',R22,'R13',R13,'R23',R23,'R33',R33);

end