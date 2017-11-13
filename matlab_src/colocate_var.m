function var = colocate_var(var,node,sur)
% =========================================================================
% colocale variable on nodes
% 
% CALL:  colocate_var(var,node,sur)
% INPUT: var - INOUT variable
%        node - 'w' colocate on w-nodes
%               'uvp' colocate on uvp-nodes
%        sur - (OPTIONAL) surface value of variable (used only for w-nodes)
%                         if not presents sur=0;
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
n = size(var);

vartmp = var;
var = zeros(n);
switch node
    case 'w'
        if(exist('sur','var'))
            var(:,:,1,:)=sur(:,:);
        else
            % assume surface value of 0
            var(:,:,1,:)=0;
        end
        var(:,:,2:n(3),:)=0.5*(vartmp(:,:,1:n(3)-1,:)+vartmp(:,:,2:n(3),:));
    case 'uvp'
        var(:,:,1:n(3)-1,:)=0.5*(vartmp(:,:,1:n(3)-1,:)+vartmp(:,:,2:n(3),:));
        var(:,:,n(3),:)=0;
end
end