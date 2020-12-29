function [] = plot4pack(PrVCR,PrJVCR,PrGCR,PrJGCR,varargin)

%==========================================================================
% CREATE FIGURE WINDOW
%-----------------------------------
fh1=figure('Units','pixels','Position',[20 40 1400 750],'Color','w');
ax1=axes('Units','pixels','Position',[100 440 550 290],'Color','none');
ax2=axes('Units','pixels','Position',[780 440 550 290],'Color','none');
ax3=axes('Units','pixels','Position',[100 75 550 290],'Color','none');
ax4=axes('Units','pixels','Position',[780 75 550 290],'Color','none');
%---AXFORMAT()
ARGS.XL = 'Number of variants considered'; 
ARGS.YL = 'P( var )';
AXFORMAT(ax1,ARGS);
ARGS.XL = 'Number of variants considered'; 
ARGS.YL = 'P( var , var )';
AXFORMAT(ax2,ARGS);
ARGS.XL = 'Number of genes considered'; 
ARGS.YL = 'P( gene )';
AXFORMAT(ax3,ARGS);
ARGS.XL = 'Number of genes considered'; 
ARGS.YL = 'P( gene , gene )';
AXFORMAT(ax4,ARGS);

fh1.Renderer = 'painters'; hold on;
%---



n = 1;
CVEC = lines(n);
ph1 = scatter(ax1,1:numel(PrVCR),  PrVCR,  500, CVEC,'Marker','.','MarkerEdgeAlpha',.4); hold on;
ph2 = scatter(ax2,1:numel(PrJVCR), PrJVCR, 500, CVEC,'Marker','.','MarkerEdgeAlpha',.4); hold on;
ph3 = scatter(ax3,1:numel(PrGCR),  PrGCR,  500, CVEC,'Marker','.','MarkerEdgeAlpha',.4); hold on;
ph4 = scatter(ax4,1:numel(PrJGCR), PrJGCR, 500, CVEC,'Marker','.','MarkerEdgeAlpha',.4); hold on;



% n = 1;
% CVEC = lines(n);
% ph1 = scatter(ax1,1:numel((PrVCR{n})), (PrVCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
% ph2 = scatter(ax2,1:numel((PrJVCR{n})), (PrJVCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
% ph3 = scatter(ax3,1:numel((PrGCR{n})), (PrGCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
% ph4 = scatter(ax4,1:numel((PrJGCR{n})), (PrJGCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;


% SAVE FIGURE TO DESKTOP
%------------------------------------------
if nargin == 6
pause(1)
FIGNAME = TBL.Properties.VariableDescriptions{'GENE'};
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/FIG_' FIGNAME '_' dt '.png']);
pause(1)
end
%------------------------------------------



end