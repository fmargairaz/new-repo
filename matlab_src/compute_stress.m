function Rij=compute_stress(ui,uj,uiuj)
% =========================================================================
% compute stress Rij=uiuj-ui*uj
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
Rij=zeros(size(ui));

for t=1:size(ui,4)
    Rij(:,:,:,t)=uiuj(:,:,:,t)-ui(:,:,:,t).*uj(:,:,:,t);
end

end