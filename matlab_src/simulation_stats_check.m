function [output] = simulation_stats_check(nx,ny,nz,lx,ly,lz,dt,z0,ustar,z_i,st_start,st_end,st_per,folder,cf,SAVEFIG_flag)
% =========================================================================
% load running average (RAV) files
% 
% CALL:  simulation_stats_check(nx,ny,nz,lx,ly,lz,dt,z0,ustar,z_i,st_start,st_end,st_per,folder,cf)
% INPUT: nx,ny,nz - grid resolution
%        lx,ly,lz - domain size
%        dt - time step
%        z0 - sruface roughness
%        ustar - friction velocity (used as normalization velocity)
%        z_i - ABL height (used as normalization lenght)
%        st_start,st_end,st_per - time step for RAV files start,end,periode
%        folder - path to data
%        cf - case name use as fieldname for the data structure
%             (will create a subfolder fif/cf for figure)
%        SAVEFIG_flag - flag to save figure (using mySave2pdf function)
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
dx=lx/nx;dy=ly/ny;dz=lz/nz; 
x=0:dx:lx-dx;y=0:dy:ly-dy;
z_w=0:dz:lz-dz;z_uvp=z_w+0.5*dz;
z_th=0:0.001:lz;u_prof_th=1/0.4*log(z_th/z0);

stat_outsteps = st_start+st_per:st_per:st_end;

data=get_data_struct_RAV(folder,nx,ny,nz,stat_outsteps);

prof=get_profile_struct(data.data3D);

Rij=compute_Rij(data.data3D);
Cij=compute_Cij(data.data3D);

runstat=load([folder,'/running_diagnostics.txt']);
t_start=runstat(1,1);
t=runstat(:,2)/ustar*z_i; 
MKE=runstat(:,4);
u_star=runstat(:,5);

avg_per=(st_start:st_end)-(t_start-1);

avg_MKE=mean(MKE(avg_per));
avg_us=mean(u_star(avg_per));

fs = {'MKE','u_star'};
cs = {avg_MKE,avg_us};
avg_run = cell2struct(cs,fs,2);

fs = {'x','y','z_w','z_uvp','t','MKE','u_star','avg_run','profiles','data2D','data3D','Rij','Cij'};
cs = {x,y,z_w,z_uvp,t,MKE,u_star,avg_run,prof,data.data2D,data.data3D,Rij,Cij};
output = cell2struct(cs,fs,2);

%==========================================================================
f1=figure('Position',[50,700,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);
plot(t,MKE)
hold all
l=line(get(gca,'Xlim'),repmat(avg_MKE,1,2));
set(l,'Color','k','LineStyle','--');
xlabel('$t [s]$')
ylabel('$MKE/u_*^2$')
xlim(t([1,end]))
h=legend('$MKE(t)$',['$\langle MKE\rangle_t=',num2str(avg_MKE),'$']);
set(h,'interpreter','latex')
set(h,'box','off','Location','SouthEast')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

%==========================================================================
f2=figure('Position',[500,700,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);
plot(t,u_star)
hold all
l=line(get(gca,'Xlim'),repmat(avg_us,1,2));
set(l,'Color','k','LineStyle','--');
xlabel('$t [s]$')
ylabel('$u_*/u_*$')
xlim(t([1,end]))
h=legend('$u_*(t)$',['$\langle u_*\rangle_t=',num2str(avg_us),'$ $u_*$']);
set(h,'interpreter','latex')
set(h,'box','off','Location','SouthEast')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

%==========================================================================
f3=figure('Position',[950,700,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);

axes(ha);
semilogx(z_uvp,prof.u,'+-')
hold on
semilogx(z_th,u_prof_th,'k--')
xlabel('$z/z_i$')
ylabel('$u$ [m/s]')
grid on
h=legend('$<u>_{xy}(z)$','$u_{Log Law}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','NorthWest')
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

%==========================================================================
f4=figure('Position',[1400,700,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);

plot(prof.tke,z_uvp,'o-')
h=legend('$<tke>_{xy}(z)$');
set(h,'interpreter','latex','Location','NorthEast')
set(h,'box','off')
xlabel('$tke/u_*^2$')
ylabel('$z/z_i$')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
set(gca,'Ylim',[0,1])

%==========================================================================
f5=figure('Position',[1400,700,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);

plot(prof.cs,z_uvp,'o-')
h=legend('$<C_S>_{xy}(z)$');
set(h,'interpreter','latex','Location','NorthEast')
set(h,'box','off')
xlabel('$C_s$')
ylabel('$z/z_i$')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
set(gca,'Ylim',[0,1])

%==========================================================================
f6=figure('Position',[1400,700,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);

plot(prof.nut,z_uvp,'o-')
h=legend('$<\nu_t>_{xy}(z)$');
set(h,'interpreter','latex','Location','NorthEast')
set(h,'box','off')
xlabel('$\nu_t$')
ylabel('$z/z_i$')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
set(gca,'Ylim',[0,1])

%==========================================================================
f7=figure('Position',[50,300,450,300]);
ha=tight_subplot(1,1,[.13 .05],[.1 .05],[.1 .05]);

axes(ha);
plot(squeeze(mean(mean(data.data3D.p,2),1)),z_uvp,'o-')
h=legend(['$<p>_{xy}(z)$']);
set(h,'interpreter','latex','Location','NorthWest')
set(h,'box','off')
xlabel('$p/u_*^2$')
ylabel('$z/z_i$')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)
set(gca,'Xlim',[0 300])
set(gca,'Ylim',[0,1])

%==========================================================================
f8=figure('Position',[50,100,450,500]);
ha=tight_subplot(1,3,[.13 .05],[.1 .1],[.1 .05]);

axes(ha(1));
plot(prof.u,z_uvp,'+-')
xlabel('$u/u_*$')
ylabel('$z/z_i$')
grid on
h=legend('$<u>_{xy}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','North')
hpos=get(h,'Position');
set(h,'Position',[hpos(1) 0.90 hpos(3) hpos(4)])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

axes(ha(2));
plot(prof.v,z_uvp,'+-')
xlabel('$v/u_*$')
grid on
h=legend('$<v>_{xy}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','North')
hpos=get(h,'Position');
set(h,'Position',[hpos(1) 0.90 hpos(3) hpos(4)])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

axes(ha(3));
plot(prof.w,z_w,'+-')
xlabel('$w/u_*$')
grid on
h=legend('$<w>_{xy}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','North')
hpos=get(h,'Position');
set(h,'Position',[hpos(1) 0.90 hpos(3) hpos(4)])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

set(ha,'Ylim',[0,1])
set(ha([2 3]),'YTickLabel','')

%==========================================================================
f9=figure('Position',[500,100,450,500]);
ha=tight_subplot(1,3,[.13 .05],[.1 .1],[.1 .05]);

axes(ha(1));
plot(prof.var_u,z_uvp,'+-')
xlabel('$(u/u_*)^2$')
ylabel('$z/z_i$')
grid on
h=legend('$<u^\prime u^\prime>_{xy}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','North')
hpos=get(h,'Position');
set(h,'Position',[hpos(1) 0.90 hpos(3) hpos(4)])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

axes(ha(2));
plot(prof.var_v,z_uvp,'+-')
xlabel('$(v/u_*)^2$')
grid on
h=legend('$<v^\prime v^\prime>_{xy}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','North')
hpos=get(h,'Position');
set(h,'Position',[hpos(1) 0.90 hpos(3) hpos(4)])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

axes(ha(3));
plot(prof.var_w,z_w,'+-')
xlabel('$(w/u_*)^2$')
grid on
h=legend('$<w^\prime w^\prime>_{xy}(z)$');
set(h,'interpreter','latex')
set(h,'box','off','Location','North')
hpos=get(h,'Position');
set(h,'Position',[hpos(1) 0.90 hpos(3) hpos(4)])
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

set(ha,'Ylim',[0,1])
set(ha([2 3]),'YTickLabel','')

%==========================================================================
f10=figure('Position',[950,100,900,500]);
plot(-prof.R13,z_w,'-')
hold all
plot(-prof.txz,z_w,'-')
plot(-Cij.C13,z_w,'-')
plot(-(prof.R13+prof.txz+Cij.C13),z_w,'-')
plot((1-z_w),z_w,'k--')
xlabel('$\tau_{xz}/u_*^2$')
ylabel('$z/z_i$')
%xlim([0 1.1])
h=legend('$<R_{xz}>_{xy}(z)$','$<\tau^{SGS}_{xz}>_{xy}(z)$',...
    '$C_{xz}(z)$','$<\tau^{SGS}_{xz}+R_{xz}>_{xy}(z)+C_{xz}(z)$',...
    'pressure forcing');
set(h,'interpreter','latex','Location','NorthEast')
set(h,'box','off')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

set(gca,'Ylim',[0,1])
set(gca,'Xlim',[0 1])

%==========================================================================
f11=figure('Position',[950,100,900,500]);
plot(-prof.R23,z_w,'+-')
hold all
plot(-prof.tyz,z_w,'+-')
plot(-Cij.C23,z_w,'+-')
plot(-(prof.R23+prof.tyz+Cij.C23),z_w,'+-')
plot(zeros(size(z_w)),z_w,'k--')
xlabel('$\tau_{yz}/u_*^2$')
ylabel('$z/z_i$')
%xlim([0 1.1])
h=legend('$<R_{yz}>_{xy}(z)$','$<\tau^{SGS}_{yz}>_{xy}(z)$',...
    '$C_{yz}(z)$','$<\tau^{SGS}_{yz}+R_{yz}>_{xy}(z)+C_{yz}(z)$',...
    'pressure forcing');
set(h,'interpreter','latex','Location','NorthEast')
set(h,'box','off')
grid on
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',2)

set(gca,'Ylim',[0,1])

%==========================================================================
if(SAVEFIG_flag) 
    system(['mkdir -p fig/',cf]);
    mySave2pdf(f1,[20 15],['fig/',cf,'/MKE_trace.pdf'])
    mySave2pdf(f2,[20 15],['fig/',cf,'/ustar_trace.pdf'])
    mySave2pdf(f3,[20 15],['fig/',cf,'/u_logprofile.pdf'])
    mySave2pdf(f4,[20 15],['fig/',cf,'/tke_profile.pdf'])
    mySave2pdf(f5,[20 15],['fig/',cf,'/cs_profile.pdf'])
    mySave2pdf(f6,[20 15],['fig/',cf,'/nut_profile.pdf'])
    mySave2pdf(f7,[20 15],['fig/',cf,'/p_profile.pdf'])
    mySave2pdf(f8,[20 15],['fig/',cf,'/vel_profile.pdf'])
    mySave2pdf(f9,[20 15],['fig/',cf,'/var_profile.pdf'])
    mySave2pdf(f10,[20 15],['fig/',cf,'/stress_XY_profile.pdf'])
    mySave2pdf(f11,[20 15],['fig/',cf,'/stress_YZ_profile.pdf'])
end
%==========================================================================
end

