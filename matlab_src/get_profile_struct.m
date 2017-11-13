function sout = get_profile_struct(s)
% =========================================================================
% compute profiles from running averages
% 
% INPUT: s - data structure 
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
u = colocate_var(s.u,'w');
v = colocate_var(s.v,'w');

var_u_colw = colocate_var(s.u2-s.u.^2,'w');
var_v_colw = colocate_var(s.v2-s.v.^2,'w');
var_w_colu = colocate_var(s.w2-s.w.^2,'uvp');

var_up = squeeze(mean(mean(s.u2-s.u.^2,2),1));
var_vp = squeeze(mean(mean(s.v2-s.v.^2,2),1));
var_wp = squeeze(mean(mean(s.w2-s.w.^2,2),1));

var_up_colw = squeeze(mean(mean(var_u_colw,2),1));
var_vp_colw = squeeze(mean(mean(var_v_colw,2),1));
var_wp_colu = squeeze(mean(mean(var_w_colu,2),1));

up = squeeze(mean(mean(s.u,2),1));
vp = squeeze(mean(mean(s.v,2),1));
wp = squeeze(mean(mean(s.w,2),1));
pp = squeeze(mean(mean(s.p,2),1));

txx = squeeze(mean(mean(s.txx,2),1));
tyy = squeeze(mean(mean(s.tyy,2),1));
tzz = squeeze(mean(mean(s.tzz,2),1));

R13 = compute_stress(u,s.w,s.uw);
St13= R13+s.txz;
R13 = squeeze(mean(mean(R13,2),1));
txz = squeeze(mean(mean(s.txz,2),1));
St13 = squeeze(mean(mean(St13,2),1));

R23 = compute_stress(v,s.w,s.vw);
St23= R23+s.tyz;
R23 = squeeze(mean(mean(R23,2),1));
tyz = squeeze(mean(mean(s.tyz,2),1));
St23 = squeeze(mean(mean(St23,2),1));

Cij=compute_Cij(s);

St13=St13+Cij.C13;
St23=St23+Cij.C23;

cs= squeeze(mean(mean(s.cs,2),1));
nut= squeeze(mean(mean(s.nut,2),1));

sout=struct('u',up,'v',vp,'w',wp,'p',pp,...
    'var_u',var_up,'var_v',var_vp,'var_w',var_wp,...
    'txx',txx,'tyy',tyy,'tzz',tzz,...
    'R13',R13,'txz',txz,'St13',St13,...
    'R23',R23,'tyz',tyz,'St23',St23,...
    'tke',0.5*(var_up+var_vp+var_wp_colu),...
    'cs',cs,'nut',nut);
end