



% Sohlh1
% In your results, you have 17 different comparisons.
% Just make sure you adjust your significant p-value by 
% dividing this value by 17


.05/17



% Pooling your data significantly reduces your statistical power

% Estrus Cycle Measurement on immune profile







WT = randn(5,17);
MU = randn(5,17);

boxplot(WT,MU)


GROUP = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

close all; boxplot(DATA',GROUP)





%%
clc; close all; clear;
%==========================================================================
% BOXPLOT USING BOXCHART()
%-----------------------------------
fh1 = figure('Units','pixels','Position',[30 40 1200 600],'Color','w');
ax1=axes('Units','normalized','Position',[.090 .13 .8 .8],'Color','none');
ax2=axes('Units','normalized','Position',[.101 .13 .8 .8],'Color','none');
%ax2 = axes('XColor','none','YColor','none','XTick',[],'YTick',[]);
axes(ax2); axis off; hold on;
%-----------------------------------


nt = 17;

WT = randn(5,nt);
MU = randn(5,nt);


axes(ax1);
ph1 = boxchart(WT, 'BoxFaceColor',[.00 .50 .50], 'BoxWidth',.15);
ylim([-4 4])
hold on;


axes(ax2)
ph2 = boxchart(MU, 'BoxFaceColor',[.50 .00 .20], 'BoxWidth',.15);
ylim([-4 4])



PVAL = zeros(nt,1);
for i = 1:nt

    [h,p,ci,stats] = ttest2(WT(:,i), MU(:,i));
    PVAL(i) = p;

end

figure; bar(sort(PVAL))



%==========================================================================
%%             BOXPLOT
%==========================================================================
clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .25 .85],'Color','w');
ax1 = axes('Position',[.09 .09 .8 .8],'Color','none');
 
    axes(ax1)
ph1=boxplot(tSNP.mD,tSNP.BRAAK,...
    'OutlierSize',1,'Symbol','.','Widths',.4,...
    'FullFactors','on','FactorGap',5,...
    'Jitter',0,'ExtremeMode','compress');
 
 
% 'Notch', 
% 	'off' | 'on' | 'marker'
% 'Labels', 
% 	{'CASE','CTRL'}
% 'Whisker', 
% 	1
% 'PlotStyle', 
% 	'traditional' | 'compact'
% 'BoxStyle', 
% 	'filled' | 'outline'
% 'Colors', 
% 	[.5 .6 .7]
% 'MedianStyle', 
% 	'line' | 'target'
% 'OutlierSize', 
% 	1
% 'Widths', 
% 	0.3
% 'ColorGroup',
% 	[.5 .6 .7;
% 	 .1 .2 .3]
% 'FactorDirection',
% 	'data' | 'list' | 'auto'
% 'FullFactors',
% 	'off' | 'on'
% 'FactorGap',
% 	5
% 'FactorSeparator',
% 	[1,2]
% 'GroupOrder',
% 	[] | string | cell
% 'DataLim', — Extreme data limits
% 	-Inf,Inf] | [-.5 .5]
% 'ExtremeMode'
% 	'clip' (default) | 'compress'
% 'Jitter'
% 	0
% 'Whisker' — Maximum whisker length
% 	1.5 (default) | positive numeric value
% 'LabelOrientation'
% 	'inline' | 'horizontal'
% 'LabelVerbosity' — Labels to display on plot
% 	'all' | 'minor' | 'majorminor'
% 'Orientation' — Plot orientation
% 	'vertical' (default) | 'horizontal'

