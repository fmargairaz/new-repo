function mySave2pdf(h,a,name,fsize)
% =========================================================================
% save figure as pdf
% 
% CALL:  mySave2pdf(h,a,name,fsize)
% INPUT: h - figure hanlde
%        a - 2D array [x,y] size of the figure in cm
%        name - file name
%        fsize - font size
%
% Fabien Margairaz, University of Utah, SLC
% =========================================================================
if(~exist('fsize'))
    fsize=14;
end

haxes = findobj(gcf, 'type', 'axes');
set(haxes,'Fontsize',fsize)

set(h,'PaperUnits','centimeters');
set(h,'PaperSize',a);
set(h,'PaperPosition', [0 0 a(1) a(2)]);

pause(1)

print(h, '-dpdf', name);

end
