function [] = plot2pack(Y1,Y2,varargin)

%==========================================================================
% CREATE FIGURE WINDOW
%-----------------------------------
fh1=figure('Units','pixels','Position',[20 90 1400 600],'Color','w');
ax1=axes('Units','pixels','Position',[100 75 550 500],'Color','none');
ax2=axes('Units','pixels','Position',[780 75 550 500],'Color','none');
%---AXFORMAT()
ARGS.XL = 'Number of genes considered'; 
ARGS.YL = 'P( gene )';
AXFORMAT(ax1,ARGS);
ARGS.XL = 'Number of genes considered'; 
ARGS.YL = 'P( gene , gene )';
AXFORMAT(ax2,ARGS);

fh1.Renderer = 'painters'; hold on;
%---



n = 1;
CVEC = lines(n);
ph1 = scatter(ax1,1:numel(Y1), Y1, 500, CVEC,'Marker','.','MarkerEdgeAlpha',.4); hold on;
ph2 = scatter(ax2,1:numel(Y2), Y2, 500, CVEC,'Marker','.','MarkerEdgeAlpha',.4); hold on;



% SAVE FIGURE TO DESKTOP
%------------------------------------------
if nargin == 3
pause(.5)
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, varargin{1});
pause(.5)
end
%------------------------------------------

end