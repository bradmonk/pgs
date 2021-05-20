%==========================================================================
%% SETUP PGS ENV
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home  = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2';
cd(P.home);
P.srip   = [P.home filesep 'PGS_SCRIPTS'];
P.json   = [P.home filesep 'PGS_JSON'];
P.gnomat = [P.home filesep 'PGS_GNOMAT'];
P.fig    = [P.home filesep 'PGS_FIGS'];
P.mat    = [P.home filesep 'PGS_MAT'];
P.csv    = [P.home filesep 'PGS_CSV'];
P.fun    = [P.home filesep 'PGS_FUNS'];
P.fip    = [P.home filesep 'PGS_FUNS/fishp'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------







%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P



PGS = load([P.mat P.f 'PGS_STEP110_OUTPUT.mat']);


GNOMAD  = PGS.GNOMAD; peek(GNOMAD)

GNOMAD = GENEijk(GNOMAD);




peek(GNOMAD); gene_report(GNOMAD)
%==========================================================================
%% LOAD TRIO DATA
%==========================================================================
clc; clearvars -except P PGS GNOMAD


%TRIO = load([P.mat '/TRIO_DATA.mat']);




%==========================================================================
%% CREATE GENE SETS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



GNO = GNOMAD;

S = GNO.PANEL_TF == 1;                 % SCREENING PANEL GENE

H = GNO.HUMAN_TARGS == 1;              % HUMAN TARGETS

M = GNO.MOUSE_TARGS == 1;              % MOUSE TARGETS



GSET.A    = GNO;                       %  ALL
GSET.S    = GNO( S , :);               %  SCREEN
GSET.H    = GNO( H , :);               %  HUMAN
GSET.M    = GNO( M , :);               %  MOUSE
GSET.HM   = GNO( (H | M) , :);         %  HUMAN | MOUSE
GSET.HMS  = GNO( (H|M|S) , :);         %  HUMAN | MOUSE | SCREEN




fprintf('\n--- VARIANT COUNT --- \n')
fprintf('UNIQUE  A    SITES:   %7.f \n',numel((GSET.A.GENE)))
fprintf('UNIQUE  S    SITES:   %7.f \n',numel((GSET.S.GENE)))
fprintf('UNIQUE  H    SITES:   %7.f \n',numel((GSET.H.GENE)))
fprintf('UNIQUE  M    SITES:   %7.f \n',numel((GSET.M.GENE)))
fprintf('UNIQUE  HM   SITES:   %7.f \n',numel((GSET.HM.GENE)))
fprintf('UNIQUE  HMS  SITES:   %7.f \n',numel((GSET.HMS.GENE)))


fprintf('\n--- GENE COUNT --- \n')
fprintf('UNIQUE  A    GENES:   %7.f \n',numel(unique(GSET.A.GENE)))
fprintf('UNIQUE  S    GENES:   %7.f \n',numel(unique(GSET.S.GENE)))
fprintf('UNIQUE  H    GENES:   %7.f \n',numel(unique(GSET.H.GENE)))
fprintf('UNIQUE  M    GENES:   %7.f \n',numel(unique(GSET.M.GENE)))
fprintf('UNIQUE  HM   GENES:   %7.f \n',numel(unique(GSET.HM.GENE)))
fprintf('UNIQUE  HMS  GENES:   %7.f \n',numel(unique(GSET.HMS.GENE)))


GSET.A.Properties.VariableDescriptions{'GENE'}    = 'A';
GSET.S.Properties.VariableDescriptions{'GENE'}    = 'S';
GSET.H.Properties.VariableDescriptions{'GENE'}    = 'H';
GSET.M.Properties.VariableDescriptions{'GENE'}    = 'M';
GSET.HM.Properties.VariableDescriptions{'GENE'}   = 'HM';
GSET.HMS.Properties.VariableDescriptions{'GENE'}  = 'HMS';





%==========================================================================
%% QUANTIFY CARRIER RATES
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET


% TBL = GSET.A;
% TBL = GSET.S;
% TBL = GSET.H;
% TBL = GSET.M;
% TBL = GSET.HM;
% TBL = GSET.HMS;


TBL = GSET.HM;





TBL = GENEijk(TBL);
TAB = TBL(TBL.GENEj == 1 ,:);


Vn = height(TBL);
Gn = height(TAB);


VTAB = table();
VTAB.VCR            = TBL.VCR;
VTAB.VCR_raw        = TBL.VCR_raw;
VTAB.VCR_popmax     = TBL.VCR_popmax;
VTAB.VCR_male       = TBL.VCR_male;
VTAB.VCR_female     = TBL.VCR_female;
VTAB.VCR_afr        = TBL.VCR_afr;
VTAB.VCR_amr        = TBL.VCR_amr;
VTAB.VCR_asj        = TBL.VCR_asj;
VTAB.VCR_eas        = TBL.VCR_eas;
VTAB.VCR_eas_jpn    = TBL.VCR_eas_jpn;
VTAB.VCR_eas_kor    = TBL.VCR_eas_kor;
VTAB.VCR_eas_oea    = TBL.VCR_eas_oea;
VTAB.VCR_fin        = TBL.VCR_fin;
VTAB.VCR_nfe        = TBL.VCR_nfe;
VTAB.VCR_nfe_bgr    = TBL.VCR_nfe_bgr;
VTAB.VCR_nfe_est    = TBL.VCR_nfe_est;
VTAB.VCR_nfe_nwe    = TBL.VCR_nfe_nwe;
VTAB.VCR_nfe_onf    = TBL.VCR_nfe_onf;
VTAB.VCR_nfe_seu    = TBL.VCR_nfe_seu;
VTAB.VCR_nfe_swe    = TBL.VCR_nfe_swe;
VTAB.VCR_oth        = TBL.VCR_oth;
VTAB.VCR_sas        = TBL.VCR_sas;
VTAB.PrVCR          = nan(Vn,1);
VTAB.PrVCR_raw      = nan(Vn,1);
VTAB.PrVCR_popmax   = nan(Vn,1);
VTAB.PrVCR_male     = nan(Vn,1);
VTAB.PrVCR_female   = nan(Vn,1);
VTAB.PrVCR_afr      = nan(Vn,1);
VTAB.PrVCR_amr      = nan(Vn,1);
VTAB.PrVCR_asj      = nan(Vn,1);
VTAB.PrVCR_eas      = nan(Vn,1);
VTAB.PrVCR_eas_jpn  = nan(Vn,1);
VTAB.PrVCR_eas_kor  = nan(Vn,1);
VTAB.PrVCR_eas_oea  = nan(Vn,1);
VTAB.PrVCR_fin      = nan(Vn,1);
VTAB.PrVCR_nfe      = nan(Vn,1);
VTAB.PrVCR_nfe_bgr  = nan(Vn,1);
VTAB.PrVCR_nfe_est  = nan(Vn,1);
VTAB.PrVCR_nfe_nwe  = nan(Vn,1);
VTAB.PrVCR_nfe_onf  = nan(Vn,1);
VTAB.PrVCR_nfe_seu  = nan(Vn,1);
VTAB.PrVCR_nfe_swe  = nan(Vn,1);
VTAB.PrVCR_oth      = nan(Vn,1);
VTAB.PrVCR_sas      = nan(Vn,1);
VTAB.PrJVCR          = nan(Vn,1);
VTAB.PrJVCR_raw      = nan(Vn,1);
VTAB.PrJVCR_popmax   = nan(Vn,1);
VTAB.PrJVCR_male     = nan(Vn,1);
VTAB.PrJVCR_female   = nan(Vn,1);
VTAB.PrJVCR_afr      = nan(Vn,1);
VTAB.PrJVCR_amr      = nan(Vn,1);
VTAB.PrJVCR_asj      = nan(Vn,1);
VTAB.PrJVCR_eas      = nan(Vn,1);
VTAB.PrJVCR_eas_jpn  = nan(Vn,1);
VTAB.PrJVCR_eas_kor  = nan(Vn,1);
VTAB.PrJVCR_eas_oea  = nan(Vn,1);
VTAB.PrJVCR_fin      = nan(Vn,1);
VTAB.PrJVCR_nfe      = nan(Vn,1);
VTAB.PrJVCR_nfe_bgr  = nan(Vn,1);
VTAB.PrJVCR_nfe_est  = nan(Vn,1);
VTAB.PrJVCR_nfe_nwe  = nan(Vn,1);
VTAB.PrJVCR_nfe_onf  = nan(Vn,1);
VTAB.PrJVCR_nfe_seu  = nan(Vn,1);
VTAB.PrJVCR_nfe_swe  = nan(Vn,1);
VTAB.PrJVCR_oth      = nan(Vn,1);
VTAB.PrJVCR_sas      = nan(Vn,1);


GTAB = table();
GTAB.GCR            = TAB.GCR;
GTAB.GCR_raw        = TAB.GCR_raw;
GTAB.GCR_popmax     = TAB.GCR_popmax;
GTAB.GCR_male       = TAB.GCR_male;
GTAB.GCR_female     = TAB.GCR_female;
GTAB.GCR_afr        = TAB.GCR_afr;
GTAB.GCR_amr        = TAB.GCR_amr;
GTAB.GCR_asj        = TAB.GCR_asj;
GTAB.GCR_eas        = TAB.GCR_eas;
GTAB.GCR_eas_jpn    = TAB.GCR_eas_jpn;
GTAB.GCR_eas_kor    = TAB.GCR_eas_kor;
GTAB.GCR_eas_oea    = TAB.GCR_eas_oea;
GTAB.GCR_fin        = TAB.GCR_fin;
GTAB.GCR_nfe        = TAB.GCR_nfe;
GTAB.GCR_nfe_bgr    = TAB.GCR_nfe_bgr;
GTAB.GCR_nfe_est    = TAB.GCR_nfe_est;
GTAB.GCR_nfe_nwe    = TAB.GCR_nfe_nwe;
GTAB.GCR_nfe_onf    = TAB.GCR_nfe_onf;
GTAB.GCR_nfe_seu    = TAB.GCR_nfe_seu;
GTAB.GCR_nfe_swe    = TAB.GCR_nfe_swe;
GTAB.GCR_oth        = TAB.GCR_oth;
GTAB.GCR_sas        = TAB.GCR_sas;
GTAB.PrGCR          = nan(Gn,1);
GTAB.PrGCR_raw      = nan(Gn,1);
GTAB.PrGCR_popmax   = nan(Gn,1);
GTAB.PrGCR_male     = nan(Gn,1);
GTAB.PrGCR_female   = nan(Gn,1);
GTAB.PrGCR_afr      = nan(Gn,1);
GTAB.PrGCR_amr      = nan(Gn,1);
GTAB.PrGCR_asj      = nan(Gn,1);
GTAB.PrGCR_eas      = nan(Gn,1);
GTAB.PrGCR_eas_jpn  = nan(Gn,1);
GTAB.PrGCR_eas_kor  = nan(Gn,1);
GTAB.PrGCR_eas_oea  = nan(Gn,1);
GTAB.PrGCR_fin      = nan(Gn,1);
GTAB.PrGCR_nfe      = nan(Gn,1);
GTAB.PrGCR_nfe_bgr  = nan(Gn,1);
GTAB.PrGCR_nfe_est  = nan(Gn,1);
GTAB.PrGCR_nfe_nwe  = nan(Gn,1);
GTAB.PrGCR_nfe_onf  = nan(Gn,1);
GTAB.PrGCR_nfe_seu  = nan(Gn,1);
GTAB.PrGCR_nfe_swe  = nan(Gn,1);
GTAB.PrGCR_oth      = nan(Gn,1);
GTAB.PrGCR_sas      = nan(Gn,1);
GTAB.PrJGCR          = nan(Gn,1);
GTAB.PrJGCR_raw      = nan(Gn,1);
GTAB.PrJGCR_popmax   = nan(Gn,1);
GTAB.PrJGCR_male     = nan(Gn,1);
GTAB.PrJGCR_female   = nan(Gn,1);
GTAB.PrJGCR_afr      = nan(Gn,1);
GTAB.PrJGCR_amr      = nan(Gn,1);
GTAB.PrJGCR_asj      = nan(Gn,1);
GTAB.PrJGCR_eas      = nan(Gn,1);
GTAB.PrJGCR_eas_jpn  = nan(Gn,1);
GTAB.PrJGCR_eas_kor  = nan(Gn,1);
GTAB.PrJGCR_eas_oea  = nan(Gn,1);
GTAB.PrJGCR_fin      = nan(Gn,1);
GTAB.PrJGCR_nfe      = nan(Gn,1);
GTAB.PrJGCR_nfe_bgr  = nan(Gn,1);
GTAB.PrJGCR_nfe_est  = nan(Gn,1);
GTAB.PrJGCR_nfe_nwe  = nan(Gn,1);
GTAB.PrJGCR_nfe_onf  = nan(Gn,1);
GTAB.PrJGCR_nfe_seu  = nan(Gn,1);
GTAB.PrJGCR_nfe_swe  = nan(Gn,1);
GTAB.PrJGCR_oth      = nan(Gn,1);
GTAB.PrJGCR_sas      = nan(Gn,1);




% REMOVE OUTLIER VCR & GCR SITES
%------------------------------------------------------------
% VCR = VTAB.VCR;
% GCR = GTAB.GCR;
% 
% close all; histogram(VCR,100)
% [VCRn,VCRe,VCRb] = histcounts(VCR,0:.01:1);
% VCRh = [VCRn' VCRe(2:end)'];
% 
% close all; histogram(GCR,100)
% [GCRn,GCRe,GCRb] = histcounts(GCR,0:.01:1);
% GCRh = [GCRn' GCRe(2:end)'];


VTAB( VTAB.VCR > .15 ,:) = [];
GTAB( GTAB.GCR > .15 ,:) = [];

%------------------------------------------------------------



VTAB_VARS = VTAB.Properties.VariableNames;
GTAB_VARS = GTAB.Properties.VariableNames;



j = 23;
k = 45;
for i = 1:22


    VTAB.(VTAB_VARS{j}) = ethnic_pr_vcr(    VTAB.(VTAB_VARS{i})  );
    VTAB.(VTAB_VARS{k}) = ethnic_pr_jvcr(   VTAB.(VTAB_VARS{i})  );


j = j+1;
k = k+1;
end






j = 23;
k = 45;
for i = 1:22


    GTAB.(GTAB_VARS{j}) = ethnic_pr_gcr(    GTAB.(GTAB_VARS{i})  );
    GTAB.(GTAB_VARS{k}) = ethnic_pr_jgcr(   GTAB.(GTAB_VARS{i})  );


j = j+1;
k = k+1;
end





%==========================================================================
%% PREP DATA FOR PRISM EXPORT
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB

%clc; disp(char(string(VTAB.Properties.VariableNames')));
%clc; disp(char(string(GTAB.Properties.VariableNames')));



VTAB_VARS = VTAB.Properties.VariableNames;
GTAB_VARS = GTAB.Properties.VariableNames;



MAX_PrVCR  = zeros(22,1);
MAX_PrJVCR = zeros(22,1);
MAX_PrGCR  = zeros(22,1);
MAX_PrJGCR = zeros(22,1);



% GET MAX CARRIER RATES
%----------------------------------------------------------------------
for i = 1:22
j = i+22;
k = i+44;


    MAX_PrVCR(i,1)   = max(VTAB.(VTAB_VARS{j}));
    MAX_PrJVCR(i,1)  = max(VTAB.(VTAB_VARS{k}));
    MAX_PrGCR(i,1)   = max(GTAB.(GTAB_VARS{j}));
    MAX_PrJGCR(i,1)  = max(GTAB.(GTAB_VARS{k}));


end






% GET MOVING AVERAGE CARRIER RATES
%----------------------------------------------------------------------

WINDOW_SIZE = 26;

for i = 1:22
j = i+22;
k = i+44;

    MMU_PrVCR(:,i)   = movmean( VTAB.(VTAB_VARS{j}), WINDOW_SIZE);
    MMU_PrJVCR(:,i)  = movmean( VTAB.(VTAB_VARS{k}), WINDOW_SIZE);
    MMU_PrGCR(:,i)   = movmean( GTAB.(GTAB_VARS{j}), WINDOW_SIZE);
    MMU_PrJGCR(:,i)  = movmean( GTAB.(GTAB_VARS{k}), WINDOW_SIZE);

end




MMU_PrVCR(1,:)   = 0;
MMU_PrJVCR(1,:)  = 0;
MMU_PrGCR(1,:)   = 0;
MMU_PrJGCR(1,:)  = 0;

MMU_PrVCR(2,:)   = .011;
MMU_PrJVCR(2,:)  = .0005;
MMU_PrGCR(2,:)   = .015;
MMU_PrJGCR(2,:)  = .001;



WINDOW_SIZE = 5;
for i = 1:22
    MMU_PrVCR(:,i)   = movmean( MMU_PrVCR(:,i), WINDOW_SIZE);
    MMU_PrJVCR(:,i)  = movmean( MMU_PrJVCR(:,i), WINDOW_SIZE);
    MMU_PrGCR(:,i)   = movmean( MMU_PrGCR(:,i), WINDOW_SIZE);
    MMU_PrJGCR(:,i)  = movmean( MMU_PrJGCR(:,i), WINDOW_SIZE);
end






LABELS = [
    "all", "raw", "popmax", "male", "female", ...
    "afr", "amr", "asj", "eas", "eas jpn", "eas kor", "eas oea", "fin", ...
    "nfe", "nfe bgr", "nfe est", "nfe nwe", "nfe onf", "nfe seu", "nfe swe",...
    "oth", "sas"];

%-----------------------------------
% REMOVE: raw, popmax, male, female
%-----------------------------------
%  1 all
% *2 raw
% *3 popmax
% *4 male
% *5 female
%  6 afr
%  7 amr
%  8 asj
%  9 eas
% 10 eas_jpn
% 11 eas_kor
% 12 eas_oea
% 13 fin
% 14 nfe
% 15 nfe_bgr
% 16 nfe_est
% 17 nfe_nwe
% 18 nfe_onf
% 19 nfe_seu
% 20 nfe_swe
% 21 oth
% 22 sas
%-----------------------------------

LABELS = LABELS([1,6:end]);

MAX_SNV_CCR = MAX_PrVCR(:,[1,6:end]);     % MAX SNV  CCR
MAX_SNV_ACR = MAX_PrJVCR(:,[1,6:end]);    % MAX SNV  ACR
MAX_GEN_CCR = MAX_PrGCR(:,[1,6:end]);     % MAX GENE CCR
MAX_GEN_ACR = MAX_PrJGCR(:,[1,6:end]);    % MAX GENE ACR


MMU_SNV_CCR = MMU_PrVCR(:,[1,6:end]);     % SMOOTHED SNV  CCR
MMU_SNV_ACR = MMU_PrJVCR(:,[1,6:end]);    % SMOOTHED SNV  ACR
MMU_GEN_CCR = MMU_PrGCR(:,[1,6:end]);     % SMOOTHED GENE CCR
MMU_GEN_ACR = MMU_PrJGCR(:,[1,6:end]);    % SMOOTHED GENE ACR



%==========================================================================
%X = MMU_SNV_CCR;
%X = MMU_SNV_ACR;
X = MMU_GEN_CCR;
%X = MMU_GEN_ACR;

%--------------------------------------------------------------------------
% FIGURE_DESCRIPTION
%--------------------------------------------------------------------------
close all;
fh1 = figure('Units','pixels','Position',[30 40 1050 700],'Color','w');
ax1 = axes('Position',[.12 .14 .82 .8],'Color','none');
xlabel('genes included'); ylabel('CCR');
ARGS.XP = -.08; ARGS.YP = -.09; AXFORMAT(ax1,ARGS);
fh1.Renderer = 'painters'; hold on;
%-----------------------------------

plot(X)
legend(LABELS,'Location','northwest','NumColumns',2)







clearvars -except P PGS GNOMAD GSET VTAB GTAB ...
    MAX_SNV_CCR MAX_SNV_ACR MAX_GEN_CCR MAX_GEN_ACR ...
    MMU_SNV_CCR MMU_SNV_ACR MMU_GEN_CCR MMU_GEN_ACR ...
    LABELS





%==========================================================================
%% PERFORM STATISTICAL ANALYSIS
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB
%addpath([P.fx '/fishp'])


HMTAB = GSET.HM;

disp(char(string(HMTAB.Properties.VariableNames')));


VARNAMES.AN = [...
"AN";
"AN_afr";
"AN_amr";
"AN_asj";
"AN_fin";
"AN_oth";
"AN_sas";
"AN_eas";
"AN_eas_jpn";
"AN_eas_kor";
"AN_eas_oea";
"AN_nfe";
"AN_nfe_bgr";
"AN_nfe_est";
"AN_nfe_nwe";
"AN_nfe_onf";
"AN_nfe_seu";
"AN_nfe_swe"];

VARNAMES.AC = [...
"AC";
"AC_afr";
"AC_amr";
"AC_asj";
"AC_fin";
"AC_oth";
"AC_sas"
"AC_eas";
"AC_eas_jpn";
"AC_eas_kor";
"AC_eas_oea";
"AC_nfe";
"AC_nfe_bgr";
"AC_nfe_est";
"AC_nfe_nwe";
"AC_nfe_onf";
"AC_nfe_seu";
"AC_nfe_swe"];

VARNAMES.HC = [...
"nhomalt";
"nhomalt_afr";
"nhomalt_amr";
"nhomalt_asj";
"nhomalt_fin";
"nhomalt_oth";
"nhomalt_sas"
"nhomalt_eas";
"nhomalt_eas_jpn";
"nhomalt_eas_kor";
"nhomalt_eas_oea";
"nhomalt_nfe";
"nhomalt_nfe_bgr";
"nhomalt_nfe_est";
"nhomalt_nfe_nwe";
"nhomalt_nfe_onf";
"nhomalt_nfe_seu";
"nhomalt_nfe_swe"];


VARNAMES.ETHNICITY = [...
"all";
"afr";
"amr";
"asj";
"fin";
"oth";
"sas";
"eas";
"eas_jpn";
"eas_kor";
"eas_oea";
"nfe";
"nfe_bgr";
"nfe_est";
"nfe_nwe";
"nfe_onf";
"nfe_seu";
"nfe_swe"];



FISHER_PVAL = array2table(zeros(size(VARNAMES.ETHNICITY,1)));
FISHER_PVAL.Properties.VariableNames = VARNAMES.ETHNICITY;
FISHER_PVAL.Properties.RowNames      = VARNAMES.ETHNICITY;

FISHER_ODDS = array2table(zeros(size(VARNAMES.ETHNICITY,1)));
FISHER_ODDS.Properties.VariableNames = VARNAMES.ETHNICITY;
FISHER_ODDS.Properties.RowNames      = VARNAMES.ETHNICITY;




for i=1:numel(VARNAMES.AN)
for j=1:numel(VARNAMES.AN)


    G1_AN   = round(mean(HMTAB.(VARNAMES.AN{i})));
    G2_AN   = round(mean(HMTAB.(VARNAMES.AN{j})));
    G1_ALTS = round(mean(HMTAB.(VARNAMES.AC{i})));
    G2_ALTS = round(mean(HMTAB.(VARNAMES.AC{j})));
    G1_REFS = G1_AN - G1_ALTS;
    G2_REFS = G2_AN - G2_ALTS;

    [PVAL, ODDS] = fishp(G1_REFS, G1_ALTS, G2_REFS, G2_ALTS);

    FISHER_PVAL.(VARNAMES.ETHNICITY{i})(j) = PVAL;
    FISHER_ODDS.(VARNAMES.ETHNICITY{i})(j) = ODDS;


end
end



close all;

ethnicity = regexprep(VARNAMES.ETHNICITY,'.+_','');

pvalues = table2array(FISHER_PVAL);

ColorDisplayData = pvalues;
ColorData        = pvalues;
ColorData(ColorData < .05) = NaN;
ColorData(ColorData > .95) = .95;
ColorData = tril(ColorData);
ColorData(ColorData == 0) = 1;

%cmap = flipud(bone);
%h = heatmap(ethnicity, ethnicity, ColorData, 'Colormap',cmap);
h = heatmap(ethnicity, ethnicity, ColorData);

cmap = h.Colormap;
cmap(end,:) = 0;
h.Colormap = cmap;

h.MissingDataLabel = 'p < .05';
h.MissingDataColor = [.9 .5 .5];

h.Title = 'P-values (Fishers Exact Test)';
% h.XLabel = 'Sizes';
% h.YLabel = 'Colors';









%% EOF