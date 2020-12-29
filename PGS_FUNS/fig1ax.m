function [fh1,ax1] = fig1ax(XLAB,YLAB)



fh1 = figure('Units','pixels','Position',[30 40 1050 700],'Color','w');

ax1 = axes('Position',[.12 .14 .82 .8],'Color','none');

xlabel(XLAB); 

ylabel(YLAB);

ARGS.XP = -.08; 

ARGS.YP = -.09; 

AXFORMAT(ax1,ARGS);

fh1.Renderer = 'painters';




end