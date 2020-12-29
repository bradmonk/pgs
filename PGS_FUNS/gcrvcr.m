function [varargout] = gcrvcr(TBL)

%==========================================================================
%% ASSIGN GENEi & GENEj & GENEk  (MUST NOT REMOVE GENES AFTER THIS FUNC)
%==========================================================================


TBL = GENEijk(TBL);




%==========================================================================
%% COMPUTE VARIANT CARRIER RATE (VCR)
%==========================================================================
% 
% VCR = (AC - HOM) ./ (.5 .* AN);
% 
% 
% Where...
% AC is the total allele count for the variant, HOM is the number of
% individuals who are homozygous for the variant, and AN is the total
% number of alleles analyzed for the variant.
% 
% AN:  count of all alleles measured (REF + ALT)
% AC:  count of ALT alleles (ALT_HOM + ALT_HET)
% HOM: count of homozygous ALT people (ALT_HOM)
%--------------------------------------------------------------------------




TBL = vcr(TBL);



%==========================================================================
%% COMPUTE GENE CARRIER RATE (GCR)
%==========================================================================
% 
% GCR(g) = 1 - prod(1 - VCR(1:k));
% 
% 
% The equation above computes the gene carrier rate (GCR) for a gene 'g',
% where VCR is the variant carrier rates. The product is computed for all 
% variants 1:k of a given gene.
% 
%--------------------------------------------------------------------------



TBL = gcr(TBL);







%==========================================================================
%% WHAT IS THE PROBABILITY SOMEONE HAS A LoF VARIANT?
%==========================================================================



PrVCR = pr_vcr(TBL);



%==========================================================================
%% PROBABILITY 2 PEOPLE BOTH HAVE THE SAME LOF VARIANT
%==========================================================================


PrJVCR = pr_jvcr(TBL);




%==========================================================================
%% PROBABILITY SOMEONE CARRIES A LOF GENE
%==========================================================================


PrGCR = pr_gcr(TBL);




%==========================================================================
%% PROBABILITY 2 PEOPLE BOTH CARRY THE SAME LOF GENE
%==========================================================================

PrJGCR = pr_jgcr(TBL);









%==========================================================================
%% PLOT CVCR JVCR CGCR JGCR
%==========================================================================





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

% [~,~] = fig1ax('Number of LoF sites considered',...
%                'P(cumulative VCR)');

% keyboard


axes(ax1)
ph1 = scatter(1:numel(nanmean(PrVCR,2)), nanmean(PrVCR,2), 900, [.9 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);
axes(ax2)
ph2 = scatter(1:numel(nanmean(PrJVCR,2)), nanmean(PrJVCR,2), 900, [.9 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);
axes(ax3)
ph3 = scatter(1:numel(nanmean(PrGCR,2)), nanmean(PrGCR,2), 900, [.9 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);
axes(ax4)
ph4 = scatter(1:numel(nanmean(PrJGCR,2)), nanmean(PrJGCR,2), 900, [.9 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);






%%

varargout = {PrVCR, PrJVCR, PrGCR, PrJGCR};

end