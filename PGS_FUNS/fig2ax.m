function [fh1,ax1,ax2] = fig2ax(ax1_xlab,ax1_ylab,ax2_xlab,ax2_ylab)


%==========================================================================
% FIGURE_DESCRIPTION
%-----------------------------------
fh1 = figure('Units','pixels','Position',[30 40 1200 600],'Color','w');
ax1=axes('Units','normalized','Position',[.09 .13 .41 .8],'Color','none');
ax2=axes('Units','normalized','Position',[.56 .13 .41 .8],'Color','none');
axes(ax1); xlabel(ax1_xlab); ylabel(ax1_ylab);
ARGS.XP = -.08; ARGS.YP = -.11; AXFORMAT(ax1,ARGS); hold on;
axes(ax2); xlabel(ax2_xlab); ylabel(ax2_ylab);
ARGS.XP = -.08; ARGS.YP = -.11; AXFORMAT(ax2,ARGS); hold on;
fh1.Renderer = 'painters'; 
%-----------------------------------


end