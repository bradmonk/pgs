%==========================================================================
%% SETUP PGS ENV
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home  = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2';
cd(P.home);
P.srip   = [P.home filesep 'PGS_SCRIPTS'];
P.json   = [P.home filesep 'PGS_JSON'];
P.gnomat = [P.home filesep 'PGS_GNOMAT'];
P.mat    = [P.home filesep 'PGS_MAT'];
P.csv    = [P.home filesep 'PGS_CSV'];
P.funs   = [P.home filesep 'PGS_FUNS'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;



%==========================================================================
%% ABOUT DAWES LEK & COOPER (2019) SETS
%==========================================================================
%{
From the article...


DAWES LEK & COOPER (2019)

Gene discovery informatics toolkit deﬁnes candidate genes for unexplained
infertility and prenatal or infantile mortality



To advance gene discovery in prenatal-lethal disorders, Dawes et al. 
created an easy-to-mine database integrating known human phenotypes with 
inheritance pattern, scores of genetic constraint, and murine and cellular 
knockout phenotypes.





Loss of Function (LoF)

Predicted LoF (pLoF) 

pLoF Intolerance (pLI)


Constraint = high pLI









To further refine our analysis, search the variants described in Clinvar of
our genes with an OMIM id. Are you supposed to find out the stats for each
variant instead of the gene as a whole? I am unsure 


CADD score: "The Combined Annotation Dependent Depletion (CADD) tool scores
the predicted deleteriousness of single nucleotide variants and
insertion/deletions variants in the human genome by integrating multiple
annotations including conservation and functional information into one
metric."
https://uswest.ensembl.org/info/genome/variation/prediction/protein_function.html
Should we look into the genes without omim annotation or for all of them?


Find out the probability of carrying a mutation in more than one gene and
of finding a partner with a mutation in any of those same genes 


We are doing our probability studies based on the OMIM database. In order
to prove that our prediction is right, we could look for those genes in the
families/trios we have sequenced and see if we see mutations in those
candidate genes and in which frequency. I just shared these files with you.
I am trying to figure out what exactly those files are


Use the CNV database to pull out the deletions found in our 784 genes (I
think that this is in addition to assess LoF mutations since according to
Svetlana this is the most common mutation)


%}

% WHAT IS THE PROBABILITY OF SOMEONE WITH MULTIPLE VARIANTS
% CADD SCORE
% INTERESTING IDEA FOR ULTIMATELY 
% COPY NUMBER VARIANT DATABASE


%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P



PGS = load([P.mat P.f 'PGS_STEP058_OUTPUT.mat']);


GNOMAD  = PGS.GNOMAD; peek(GNOMAD)





%==========================================================================
%% EXAMINE GNOMAD DATASET
%==========================================================================
% 
% GNOMAD DATASET COLUMN REFERENCE
%----------------------------------------------------------------
% GNOMAD COLUMNS    T.GENEi               -    T.isSV
% PLI    COLUMNS    T.pLI                 -    T.oe_lof
% DAWES  COLUMNS    T.hgnc_id             -    T.isdawes
% MK1    COLUMNS    T.MKO1_TARGETS        -    T.placenta
% MK2    COLUMNS    T.HOMOLO_GENE_ID      -    T.MKO2_TARGETS
% SCR    COLUMNS    T.PANEL               -    T.PANEL_TF
% TARG   COLUMNS    T.HAS_OMIM            -    T.MOUSE2_TARGETS
% OMIM   COLUMNS    T.MIMNumber           -    T.MIMNumber
% CLIN   COLUMNS    T.CLINVAR_PATHOGENIC  -    T.CLINVAR_BENIGN
% CADD   COLUMNS    T.CADD_RAW            -    T.CADD_PHRED 
% (297 COLUMNS)
%----------------------------------------------------------------
clc; clearvars -except P PGS GNOMAD





sum(GNOMAD.PANEL_TF) / height(GNOMAD)
sum(GNOMAD.HAS_OMIM) / height(GNOMAD)
sum(~isnan(GNOMAD.MIMNumber)) / height(GNOMAD)
sum(~isnan(GNOMAD.CLINVAR_PATHOGENIC)) / height(GNOMAD)
sum(~isnan(GNOMAD.CADD_RAW)) / height(GNOMAD)









%==========================================================================
%% CREATE GENE SETS
%==========================================================================
clc; clearvars -except P PGS GNOMAD


p = GNOMAD.PANEL_TF == 1;          % SCREEN PANEL   TARGETS (PAN) (P)
d = GNOMAD.DAWES_TARGETS == 1;     % DAWES NON-OMIM TARGETS (DAW) (D)
h = GNOMAD.HUMAN_TARGETS == 1;     % DAWES HUMAN    TARGETS (HUM) (H)
m = GNOMAD.MOUSE1_TARGETS == 1;    % MOUSE LIST-1   TARGETS (MK1) (M)
k = GNOMAD.MOUSE2_TARGETS == 1;    % MOUSE LIST-2   TARGETS (MK2) (K)



GSET.GNO    = GNOMAD;                           %  GNO
GSET.PAN    = GNOMAD( p , :);                   %  PAN
GSET.nPAN   = GNOMAD( ~p , :);                  % ~PAN
GSET.DAW    = GNOMAD( d , :);                   %  DAW
GSET.nDAW   = GNOMAD( ~d , :);                  % ~DAW
GSET.HUM    = GNOMAD( h , :);                   %  HUM
GSET.nHUM   = GNOMAD( ~h , :);                  % ~HUM
GSET.MK1    = GNOMAD( m , :);                   %  MK1
GSET.nMK1   = GNOMAD( ~m , :);                  % ~MK1
GSET.MK2    = GNOMAD( k , :);                   %  MK2
GSET.nMK2   = GNOMAD( ~k , :);                  % ~MK2

GSET.DM     = GNOMAD( d | m , :);               %  D | M
GSET.DK     = GNOMAD( d | k , :);               %  D | K
GSET.DnP    = GNOMAD( d & ~p , :);              %  D & ~P
GSET.DMKnP  = GNOMAD( (d|m|k) & ~p , :);        % (D | M | K) & ~P
GSET.DMK    = GNOMAD( (d|m|k) , :);             % (D | M | K)

GSET.HM     = GNOMAD( h | m , :);               %  H | M
GSET.HK     = GNOMAD( h | k , :);           	%  H | K
GSET.HnP    = GNOMAD( h & ~p , :);            	%  H & ~P
GSET.HMKnP  = GNOMAD( (h|m|k) & ~p , :);        % (H | M | K) & ~P
GSET.HMK    = GNOMAD( (h|m|k) , :);             % (H | M | K)








%==========================================================================
%% GENE SET REPORTS
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET


fprintf('UNIQUE  GNO    GENES:   %7.f \n',numel(unique(GSET.GNO.GENE)))
fprintf('UNIQUE  PAN    GENES:   %7.f \n',numel(unique(GSET.PAN.GENE)))
fprintf('UNIQUE nPAN    GENES:   %7.f \n',numel(unique(GSET.nPAN.GENE)))
fprintf('UNIQUE  DAW    GENES:   %7.f \n',numel(unique(GSET.DAW.GENE)))
fprintf('UNIQUE nDAW    GENES:   %7.f \n',numel(unique(GSET.nDAW.GENE)))
fprintf('UNIQUE  HUM    GENES:   %7.f \n',numel(unique(GSET.HUM.GENE)))
fprintf('UNIQUE nHUM    GENES:   %7.f \n',numel(unique(GSET.nHUM.GENE)))
fprintf('UNIQUE  MK1    GENES:   %7.f \n',numel(unique(GSET.MK1.GENE)))
fprintf('UNIQUE nMK1    GENES:   %7.f \n',numel(unique(GSET.nMK1.GENE)))
fprintf('UNIQUE  MK2    GENES:   %7.f \n',numel(unique(GSET.MK2.GENE)))
fprintf('UNIQUE nMK2    GENES:   %7.f \n',numel(unique(GSET.nMK2.GENE)))

fprintf('UNIQUE  DM     GENES:   %7.f \n',numel(unique(GSET.DM.GENE)))
fprintf('UNIQUE  DK     GENES:   %7.f \n',numel(unique(GSET.DK.GENE)))
fprintf('UNIQUE  DnP    GENES:   %7.f \n',numel(unique(GSET.DnP.GENE)))
fprintf('UNIQUE  DMKnP  GENES:   %7.f \n',numel(unique(GSET.DMKnP.GENE)))
fprintf('UNIQUE  DMK    GENES:   %7.f \n',numel(unique(GSET.DMK.GENE)))

fprintf('UNIQUE  HM     GENES:   %7.f \n',numel(unique(GSET.HM.GENE)))
fprintf('UNIQUE  HK     GENES:   %7.f \n',numel(unique(GSET.HK.GENE)))
fprintf('UNIQUE  HnP    GENES:   %7.f \n',numel(unique(GSET.HnP.GENE)))
fprintf('UNIQUE  HMKnP  GENES:   %7.f \n',numel(unique(GSET.HMKnP.GENE)))
fprintf('UNIQUE  HMK    GENES:   %7.f \n',numel(unique(GSET.HMK.GENE)))





GSET.GNO.Properties.VariableDescriptions{'GENE'}    = 'GNO';
GSET.PAN.Properties.VariableDescriptions{'GENE'}    = 'PAN';
GSET.nPAN.Properties.VariableDescriptions{'GENE'}   = 'nPAN';
GSET.DAW.Properties.VariableDescriptions{'GENE'}    = 'DAW';
GSET.nDAW.Properties.VariableDescriptions{'GENE'}   = 'nDAW';
GSET.HUM.Properties.VariableDescriptions{'GENE'}    = 'HUM';
GSET.nHUM.Properties.VariableDescriptions{'GENE'}   = 'nHUM';
GSET.MK1.Properties.VariableDescriptions{'GENE'}    = 'MK1';
GSET.nMK1.Properties.VariableDescriptions{'GENE'}   = 'nMK1';
GSET.MK2.Properties.VariableDescriptions{'GENE'}    = 'MK2';
GSET.nMK2.Properties.VariableDescriptions{'GENE'}   = 'nMK2';

GSET.DM.Properties.VariableDescriptions{'GENE'}     = 'DM';
GSET.DK.Properties.VariableDescriptions{'GENE'}     = 'DK';
GSET.DnP.Properties.VariableDescriptions{'GENE'}    = 'DnP';
GSET.DMKnP.Properties.VariableDescriptions{'GENE'}  = 'DMKnP';
GSET.DMK.Properties.VariableDescriptions{'GENE'}    = 'DMK';

GSET.HM.Properties.VariableDescriptions{'GENE'}     = 'HM';
GSET.HK.Properties.VariableDescriptions{'GENE'}     = 'HK';
GSET.HnP.Properties.VariableDescriptions{'GENE'}    = 'HnP';
GSET.HMKnP.Properties.VariableDescriptions{'GENE'}  = 'HMKnP';
GSET.HMK.Properties.VariableDescriptions{'GENE'}    = 'HMK';









%==========================================================================
%% DETERMINE WHICH PANEL GENES ARE LETHAL
%==========================================================================
% 
% GNOMAD DATASET COLUMN REFERENCE
%----------------------------------------------------------------
% GNOMAD COLUMNS    T.GENEi               -    T.isSV
% PLI    COLUMNS    T.pLI                 -    T.oe_lof
% DAWES  COLUMNS    T.hgnc_id             -    T.isdawes
% MK1    COLUMNS    T.MKO1_TARGETS        -    T.placenta
% MK2    COLUMNS    T.HOMOLO_GENE_ID      -    T.MKO2_TARGETS
% SCR    COLUMNS    T.PANEL               -    T.PANEL_TF
% TARG   COLUMNS    T.HAS_OMIM            -    T.MOUSE2_TARGETS
% OMIM   COLUMNS    T.MIMNumber           -    T.MIMNumber
% CLIN   COLUMNS    T.CLINVAR_PATHOGENIC  -    T.CLINVAR_BENIGN
% CADD   COLUMNS    T.CADD_RAW            -    T.CADD_PHRED 
% (297 COLUMNS)
%----------------------------------------------------------------
% 
% TBL = GSET.GNO;
% TBL = GSET.PAN; 
% TBL = GSET.nPAN;
% TBL = GSET.DnP; 
% TBL = GSET.DMKnP;
% TBL = GSET.HnP; 
% TBL = GSET.HMKnP;
% TBL = GSET.HMK;
%---
% TBL = GSET.PAN;
% TBL = GSET.HMKnP;
% TBL = GSET.HMK;
%----------------------------------------------------------------
clc; clearvars -except P PGS GNOMAD GSET


TBL1 = GSET.PAN; 

TBL2 = TBL1(TBL1.GENEj ,:);


writetable(TBL1,'/Users/bradleymonk/Desktop/TBL1.csv')
writetable(TBL2,'/Users/bradleymonk/Desktop/TBL2.csv')




%==========================================================================
%% QUANTIFY CARRIER RATES
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET


% TBL = GSET.GNO;
% TBL = GSET.PAN; 
% TBL = GSET.nPAN;
% TBL = GSET.DAW;
% TBL = GSET.nDAW; 
% TBL = GSET.HUM;
% TBL = GSET.nHUM; 
% TBL = GSET.MK1; 
% TBL = GSET.nMK1;
% TBL = GSET.MK2; 
% TBL = GSET.nMK2;
% TBL = GSET.DM; 
% TBL = GSET.DK; 
% TBL = GSET.DnP; 
% TBL = GSET.DMKnP;
% TBL = GSET.DMK;
% TBL = GSET.HM; 
% TBL = GSET.HK; 
% TBL = GSET.HnP; 
% TBL = GSET.HMKnP;
% TBL = GSET.HMK;


% TBL = GSET.PAN;
% TBL = GSET.HMKnP;
TBL = GSET.HMK;



TBL    = GENEijk(TBL); 
TBL    = vcr(TBL); 
TBL    = gcr(TBL);
PrVCR  = nanmean(pr_vcr(TBL),2);
PrJVCR = nanmean(pr_jvcr(TBL),2);
PrGCR  = nanmean(pr_gcr(TBL),2);
PrJGCR = nanmean(pr_jgcr(TBL),2);


% PctJGCR = PrJGCR .* 100
% x = (1:numel(PrJGCR))';



close all;
plot4pack( PrVCR, PrJVCR, PrGCR, PrJGCR )

clc;
fprintf('SITES:	%9.0f \n' , height(TBL))
fprintf('GENES: %9.0f \n' , numel(unique(TBL.GENE)))
disp(' ');
fprintf('max PrVCR:  %9.4f \n' , max(PrVCR) )
fprintf('max PrJVCR: %9.4f \n' , max(PrJVCR) )
fprintf('max PrGCR:  %9.4f \n' , max(PrGCR) )
fprintf('max PrJGCR: %9.4f \n' , max(PrJGCR) )














%==========================================================================
%% QUANTIFY CARRIER RATES FOR GENES WITH OMIM NUMBER
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET


% TBL = GSET.PAN;
% TBL = GSET.HMKnP;
TBL = GSET.HMK;



i = isnan(TBL.MIMNumber);

TBL(i,:) = [];






TBL    = GENEijk(TBL); 
TBL    = vcr(TBL); 
TBL    = gcr(TBL);
PrVCR  = nanmean(pr_vcr(TBL),2);
PrJVCR = nanmean(pr_jvcr(TBL),2);
PrGCR  = nanmean(pr_gcr(TBL),2);
PrJGCR = nanmean(pr_jgcr(TBL),2);

PctJGCR = PrJGCR .* 100
x = (1:numel(PrJGCR))';



close all;
plot4pack( PrVCR, PrJVCR, PrGCR, PrJGCR )

clc;
fprintf('SITES:	%9.0f \n' , height(TBL))
fprintf('GENES: %9.0f \n' , numel(unique(TBL.GENE)))
disp(' ');
fprintf('max PrVCR:  %9.4f \n' , max(PrVCR) )
fprintf('max PrJVCR: %9.4f \n' , max(PrJVCR) )
fprintf('max PrGCR:  %9.4f \n' , max(PrGCR) )
fprintf('max PrJGCR: %9.4f \n' , max(PrJGCR) )



























%==========================================================================
%% SIMULATE OFFSPRING CARRIER RATES
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET



% TBL = GSET.PAN;
% TBL = GSET.HMKnP;
TBL = GSET.HMK;



TBL    = GENEijk(TBL); 
TBL    = vcr(TBL); 
TBL    = gcr(TBL);
PrVCR  = nanmean(pr_vcr(TBL),2);
PrJVCR = nanmean(pr_jvcr(TBL),2);
PrGCR  = nanmean(pr_gcr(TBL),2);
PrJGCR = nanmean(pr_jgcr(TBL),2);

close all;
plot4pack( PrVCR, PrJVCR, PrGCR, PrJGCR )



PrVCR = 1 - prod(1 - TBL.VCR);

PrJVCR = 1 - prod(1 - TBL.VCR.^2);

GCR = TBL.GCR(TBL.GENEj == 1).^2;
PrJGCR = 1 - prod(1 - GCR);


p1 = PrJGCR
p2 = PrJGCR * PrJGCR
p3 = PrJGCR * PrJGCR * PrJGCR
p4 = PrJGCR * PrJGCR * PrJGCR * PrJGCR
p5 = PrJGCR * PrJGCR * PrJGCR * PrJGCR * PrJGCR

ppct = [p1 p2 p3 p4 p5]' .* 100
plog = -log([p1 p2 p3 p4 p5]')


c1 = 1 - (.75)
c2 = 1 - (.75 * .75)
c3 = 1 - (.75 * .75 * .75)
c4 = 1 - (.75 * .75 * .75 * .75)
c5 = 1 - (.75 * .75 * .75 * .75 * .75)

cpct = [c1 c2 c3 c4 c5]' .* 100
clog = -log([c1 c2 c3 c4 c5]')


%==========================================================================
%% RUN CARRIER RATES BASED ON pLI
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET
%--------------------------------------------------------------------------
% Probability of loss of function intolerance (pLI), for predicted 
% loss-of-function (pLoF) variation. The higher the pLI, the higher the 
% probability that a gene is intolerant to loss-of-function.
% 
% pLI closer to 1 indicates that the gene or transcript cannot tolerate 
% protein truncating variation (nonsense, splice acceptor and splice 
% donor variation). The gnomAD team recommends transcripts with a 
% pLI >= 0.9 for the set of transcripts extremely intolerant to
% truncating variants. pLI is based on the idea that transcripts can 
% be classified into three categories:
% 
%   > null: heterozygous tolerant; homozygous tolerant
% 
%   > recessive: heterozygous tolerant; homozygous intolerant
% 
%   > haploinsufficient: heterozygous intolerant; homozygous intolerant
% 
% An expectation-maximization algorithm was then used to
% assign a probability of belonging in each class to each gene or
% transcript. pLI is the probability of belonging in the haploinsufficient
% class.
%--------------------------------------------------------------------------


%TBL = GSET.GNOMAD;
%TBL = GSET.PAN; 
%TBL = GSET.NPAN;
%TBL = GSET.MK1; 
%TBL = GSET.NMK1;
%TBL = GSET.DAW;
%TBL = GSET.NDAW; 
%TBL = GSET.DAW_NPAN; 
%TBL = GSET.MK1_DAW;
%TBL = GSET.MK1_DAW_NPAN;
TBL = GSET.MKB_DAW_NPAN;


min(TBL.pLI)
max(TBL.pLI)


t = linspace(0,1,11);
for i = 1:numel(t)
disp(i)

    TBL(TBL.pLI < t(i) ,:) = [];

    TBL = GENEijk(TBL); 
    TBL = vcr(TBL); 
    TBL = gcr(TBL);
    PrVCR{i}  = nanmean(pr_vcr(TBL),2);
    PrJVCR{i} = nanmean(pr_jvcr(TBL),2);
    PrGCR{i}  = nanmean(pr_gcr(TBL),2);
    PrJGCR{i} = nanmean(pr_jgcr(TBL),2);

end


%%


PrVCR_max  = [t' cellfun(@max,PrVCR)'];
PrJVCR_max = [t' cellfun(@max,PrJVCR)'];
PrGCR_max  = [t' cellfun(@max,PrGCR)'];
PrJGCR_max = [t' cellfun(@max,PrJGCR)'];


%close all; gcrvcr(TBL);


% SAVE FIGURE TO DESKTOP
%------------------------------------------
% pause(1)
% FIGNAME = TBL.Properties.VariableDescriptions{'GENE'};
% set(gcf, 'PaperPositionMode', 'auto');
% dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
% saveas(gcf, ['/Users/bradleymonk/Desktop/FIG_' FIGNAME '_' dt '.png']);
% pause(1)
%------------------------------------------



%==========================================================================
% CREATE FIGURE WINDOW
%-----------------------------------
close all;
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
n = numel(PrVCR);

CVEC = lines(n);
for i = 1:n
%axes(ax1)
ph1 = scatter(ax1,1:numel((PrVCR{i})), (PrVCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
%axes(ax2)
ph2 = scatter(ax2,1:numel((PrJVCR{i})), (PrJVCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
%axes(ax3)
ph3 = scatter(ax3,1:numel((PrGCR{i})), (PrGCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
%axes(ax4)
ph4 = scatter(ax4,1:numel((PrJGCR{i})), (PrJGCR{i}), 500, CVEC(i,:),'Marker','.','MarkerEdgeAlpha',.4); hold on;
end

% SAVE FIGURE TO DESKTOP
%------------------------------------------
pause(1)
FIGNAME = TBL.Properties.VariableDescriptions{'GENE'};
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/FIG_' FIGNAME '_' dt '.png']);
pause(1)
%------------------------------------------



%==========================================================================
%% RUN CARRIER RATE FUNCTIONS
%==========================================================================
clc; clearvars -except P PGS GNOMAD GAD DAW MK1 MK2 GUO GSET

%TBL = GSET.GNOMAD;
%TBL = GSET.PAN; 
%TBL = GSET.NPAN;
%TBL = GSET.MK1; 
%TBL = GSET.NMK1;
%TBL = GSET.DAW;
%TBL = GSET.NDAW; 
%TBL = GSET.MK1_DAW;



for n = 1:8


switch n
    case 1
        TBL = GSET.GNOMAD;
    case 2
        TBL = GSET.PAN;
    case 3
        TBL = GSET.NPAN;
    case 4
        TBL = GSET.MK1;
    case 5
        TBL = GSET.NMK1;
    case 6
        TBL = GSET.DAW;
    case 7
        TBL = GSET.NDAW;
    case 8
        TBL = GSET.MK1_DAW;
end



% TBL = GENEijk(TBL); 
% TBL = vcr(TBL); 
% TBL = gcr(TBL);

close all; gcrvcr(TBL);


% SAVE FIGURE TO DESKTOP
%------------------------------------------
pause(1)
FIGNAME = TBL.Properties.VariableDescriptions{'GENE'};
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/FIG_' FIGNAME '_' dt '.png']);
pause(1)
%------------------------------------------



end








%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD GAD DAW MK1 MK2 GUO GSET



save([P.mat P.f
'PGS_STEP057_OUTPUT.mat'],'GNOMAD','DAW','MK1','GUO','GSET');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP054_OUTPUT.csv'])

