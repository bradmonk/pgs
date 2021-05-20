%==========================================================================
%% SETUP PGS ENV
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home  = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2';
P.desk  = '/Users/bradleymonk/Desktop';
cd(P.home);
P.srip   = [P.home filesep 'PGS_SCRIPTS'];
P.json   = [P.home filesep 'PGS_JSON'];
P.gnomat = [P.home filesep 'PGS_GNOMAT'];
P.mat    = [P.home filesep 'PGS_MAT/run02'];
P.csv    = [P.home filesep 'PGS_CSV'];
P.fun    = [P.home filesep 'PGS_FUNS/run02'];
P.fig    = [P.home filesep 'PGS_FIGS'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------







%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P



PGS = load([P.mat P.f 'PGS_STEP059_OUTPUT.mat']);


GNOMAD  = PGS.GNOMAD; peek(GNOMAD)

GNOMAD = GENEijk(GNOMAD);





peek(GNOMAD); gene_report(GNOMAD)
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
GSET.HMP  = GNO( (H|M|S) , :);         %  HUMAN | MOUSE | SCREEN




fprintf('UNIQUE  A    GENES:   %7.f \n',numel(unique(GSET.A.GENE)))
fprintf('UNIQUE  S    GENES:   %7.f \n',numel(unique(GSET.S.GENE)))
fprintf('UNIQUE  H    GENES:   %7.f \n',numel(unique(GSET.H.GENE)))
fprintf('UNIQUE  M    GENES:   %7.f \n',numel(unique(GSET.M.GENE)))
fprintf('UNIQUE  HM   GENES:   %7.f \n',numel(unique(GSET.HM.GENE)))
fprintf('UNIQUE  HMP  GENES:   %7.f \n',numel(unique(GSET.HMP.GENE)))


GSET.A.Properties.VariableDescriptions{'GENE'}    = 'A';
GSET.S.Properties.VariableDescriptions{'GENE'}    = 'S';
GSET.H.Properties.VariableDescriptions{'GENE'}    = 'H';
GSET.M.Properties.VariableDescriptions{'GENE'}    = 'M';
GSET.HM.Properties.VariableDescriptions{'GENE'}   = 'HM';
GSET.HMP.Properties.VariableDescriptions{'GENE'}  = 'HMP';





%==========================================================================
%% QUANTIFY CARRIER RATES
%==========================================================================
clc; clearvars -except P PGS GNOMAD GSET


% TBL = GSET.A;
% TBL = GSET.S;
% TBL = GSET.H;
% TBL = GSET.M;
% TBL = GSET.HM;
% TBL = GSET.HMP;


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
%% MULTIVARIANT CARRIER RATES
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB



VCR = VTAB.VCR;


% REMOVE EXTREMELY RARE VARIANTS (VCR < .0002)
%------------------------------------------------------------
% close all; histogram(VCR,0:.00001:.0005)
% [VCRn,VCRe,VCRb] = histcounts(VCR,0:.00001:.0005);
% VCRh = [VCRn', VCRe(2:end)'];
% VCR(VCR < .0002) = [];
%------------------------------------------------------------


N_BOOTS   = 100;
N_PEOPLE  = 100000;
N_VARS    = size(VCR,1);





ODDS = round( N_PEOPLE ./ (1 ./ VCR));



VMX = zeros(N_VARS, N_PEOPLE);
BMX = zeros(N_BOOTS, N_PEOPLE);


for j = 1:N_BOOTS
disp(j)

    VMX(:) = 0;

    for i = 1:N_VARS

        VMX(i, randperm(N_PEOPLE,ODDS(i)) ) = 1;

    end

    BMX(j,:) = sum(VMX,1);

end


clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB BMX


Normalization={'count';'probability';'countdensity';'pdf';'cumcount';'cdf'};
close all; ph = histogram(BMX,-.5:1:10.5,'Normalization',Normalization{2});

AAA = ph.Values;
AAB = ph.BinEdges;







%==========================================================================
%% MULTIGENE CARRIER RATES
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB



GCR = GTAB.GCR;

N_BOOTS   = 100;
N_PEOPLE  = 100000;
N_GENES    = size(GCR,1);


ODDS = round( N_PEOPLE ./ (1 ./ GCR));


GMX = zeros(N_GENES, N_PEOPLE);
BMX = zeros(N_BOOTS, N_PEOPLE);


for j = 1:N_BOOTS
disp(j)

    GMX(:) = 0;

    for i = 1:N_GENES

        GMX(i, randperm(N_PEOPLE,ODDS(i)) ) = 1;

    end

    BMX(j,:) = sum(GMX,1);

end


clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB BMX


Normalization={'count';'probability';'countdensity';'pdf';'cumcount';'cdf'};
close all; ph = histogram(BMX,-.5:1:10.5,'Normalization',Normalization{2});

AAA = ph.Values;
AAB = ph.BinEdges;



%==========================================================================
%% MULTIGENE CARRIER RATES (ALTERNATIVE METHOD - NO LOOPS - BIG RAM)
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB


GCR = GTAB.GCR;

N_BOOTS   = 50;
N_PEOPLE  = 10000;
N_GENES    = size(GCR,1);

GCRmx = repmat(GCR,1,N_PEOPLE);

RMX = rand(N_GENES, N_PEOPLE, N_BOOTS) < GCRmx;
SMX = squeeze(sum(RMX));

Normalization={'count';'probability';'countdensity';'pdf';'cumcount';'cdf'};
close all; ph = histogram(SMX,-.5:1:10.5,'Normalization',Normalization{2});
AAA = ph.Values;
AAB = ph.BinEdges;






%==========================================================================
%% MULTIVARIANT AT-RISK COUPLE RATE
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB



VCR = VTAB.VCR;

N_BOOTS   = 20;
N_PEOPLE  = 10000;
N_VARS    = size(VCR,1);

GCRmx = repmat(VCR,1,N_PEOPLE);

RMX = rand(N_VARS, N_PEOPLE, N_BOOTS) < GCRmx;




ACR = zeros(1,N_PEOPLE);

k = 1;
for i = 1:N_BOOTS
disp(i)
for j = 1:N_BOOTS

    if i ~= j

        CMX = ( RMX(:,:,i) + RMX(:,:,j) ) > 1;
        ACR(k,:) = sum(CMX);
        k = k+1;

    end

end
end





Normalization={'count';'probability';'countdensity';'pdf';'cumcount';'cdf'};
close all; ph = histogram(ACR,.5:1:10.5,'Normalization',Normalization{2});
AAA = ph.Values;
AAB = ph.BinEdges;
fprintf('\n P = %0.20f \n', AAA .* 100)
AAC = ceil(AAA .* 10000);








%==========================================================================
%% MULTIGENE AT-RISK COUPLE RATE
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB



GCR = GTAB.GCR;

N_BOOTS   = 20;
N_PEOPLE  = 10000;
N_GENES    = size(GCR,1);

GCRmx = repmat(GCR,1,N_PEOPLE);

RMX = rand(N_GENES, N_PEOPLE, N_BOOTS) < GCRmx;




ACR = zeros(1,N_PEOPLE);

k = 1;
for i = 1:N_BOOTS
disp(i)
for j = 1:N_BOOTS

    if i ~= j

        CMX = ( RMX(:,:,i) + RMX(:,:,j) ) > 1;
        ACR(k,:) = sum(CMX);
        k = k+1;

    end

end
end





Normalization={'count';'probability';'countdensity';'pdf';'cumcount';'cdf'};
close all; ph = histogram(ACR,.5:1:10.5,'Normalization',Normalization{2});
AAA = ph.Values;
AAB = ph.BinEdges;
fprintf('\n P = %0.20f \n', AAA .* 100)
AAC = ceil(AAA .* 100000);




%==========================================================================
%% PLOT CARRIER RATES
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB


VTAB_VARS = VTAB.Properties.VariableNames;
GTAB_VARS = GTAB.Properties.VariableNames;


for i = 1:22
j = i+22;
k = i+44;
close all;


    PrVCR   = VTAB.(VTAB_VARS{j});
    PrJVCR  = VTAB.(VTAB_VARS{k});
    PrGCR   = GTAB.(GTAB_VARS{j});
    PrJGCR  = GTAB.(GTAB_VARS{k});



    
    
    plot2pack(PrGCR, PrJGCR,  [P.fig '/PGS_MAF50_noOL_' GTAB_VARS{k} '.png'])
    %plot4pack( PrVCR, PrJVCR, PrGCR, PrJGCR )


    fprintf('\n\nVARIABLES: %s, %s, %s, %s \n\n',...
        VTAB_VARS{j}, VTAB_VARS{k}, GTAB_VARS{j}, GTAB_VARS{k})
    fprintf('SITES:	%9.0f \n' , height(VTAB))
    fprintf('GENES: %9.0f \n' , height(GTAB))
    disp(' ');
    fprintf('max PrVCR:  %9.4f \n' , max(PrVCR) )
    fprintf('max PrJVCR: %9.4f \n' , max(PrJVCR) )
    fprintf('max PrGCR:  %9.4f \n' , max(PrGCR) )
    fprintf('max PrJGCR: %9.4f \n' , max(PrJGCR) )


end






%==========================================================================
%% PREP DATA FOR PRISM EXPORT
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB


VTAB_VARS = VTAB.Properties.VariableNames;
GTAB_VARS = GTAB.Properties.VariableNames;

MAX_PrVCR  = zeros(22,1);
MAX_PrJVCR = zeros(22,1);
MAX_PrGCR  = zeros(22,1);
MAX_PrJGCR = zeros(22,1);

for i = 1:22
j = i+22;
k = i+44;


    MAX_PrVCR(i,1)   = max(VTAB.(VTAB_VARS{j}));
    MAX_PrJVCR(i,1)  = max(VTAB.(VTAB_VARS{k}));
    MAX_PrGCR(i,1)   = max(GTAB.(GTAB_VARS{j}));
    MAX_PrJGCR(i,1)  = max(GTAB.(GTAB_VARS{k}));


end




for i = 1:22
j = i+22;
k = i+44;


    MOVMU_PrVCR(:,i)   = movmean(VTAB.(VTAB_VARS{j}),26);
    MOVMU_PrJVCR(:,i)  = movmean(VTAB.(VTAB_VARS{k}),26);
    MOVMU_PrGCR(:,i)   = movmean(GTAB.(GTAB_VARS{j}),26);
    MOVMU_PrJGCR(:,i)  = movmean(GTAB.(GTAB_VARS{k}),26);


end



MM_PrVCR  = MOVMU_PrVCR(1:13:end,:);
MM_PrJVCR = MOVMU_PrJVCR(1:13:end,:);
MM_PrGCR  = MOVMU_PrGCR(1:13:end,:);
MM_PrJGCR = MOVMU_PrJGCR(1:13:end,:);



ASIAN_PrVCR  = MOVMU_PrVCR(1:26:end,:);
ASIAN_PrJVCR = MOVMU_PrJVCR(1:26:end,:);
ASIAN_PrGCR  = MOVMU_PrGCR(1:26:end,:);
ASIAN_PrJGCR = MOVMU_PrJGCR(1:26:end,:);




%{
MAX_PrGCR = MAX_PrGCR';
MAX_PrJGCR = MAX_PrJGCR';

disp(char(string(GTAB_VARS')))


 1 GCR           
 2 GCR raw       
 3 GCR popmax    
 4 GCR male      
 5 GCR female    
 6 AFR
 7 AMR
 8 ASJ
 9 EAS
10 EAS JPN
11 EAS KOR
12 EAS OEA
13 FIN
14 NFE
15 NFE BGR
16 NFE EST
17 NFE NWE
18 NFE ONF
19 NFE SEU
20 NFE SWE
21 OTH 
22 SAS

AAA = ["afr", "amr", "asj", "eas", "eas jpn", "eas kor", "eas oea", "fin", ...
       "nfe", "nfe bgr", "nfe est", "nfe nwe", "nfe onf", "nfe seu", "nfe swe",...
       "oth", "sas"]


AAB = (1:676)';

AAB = (0:26:676)';

AAC = zeros(1,17)
%}









% PrVCR  = ethnic_pr_vcr(   TBL.VCR  );
% 
% PrJVCR = ethnic_pr_jvcr(  TBL.VCR  );
% 
% PrGCR  = ethnic_pr_gcr(   TBL.GCR(TBL.GENEj == 1)  );
% 
% PrJGCR = ethnic_pr_jgcr(  TBL.GCR(TBL.GENEj == 1)  );







%==========================================================================
%% PLOT VCR & GCR DISTRIBUTIONS
%==========================================================================
clc; close all; clearvars -except P PGS GNOMAD GSET VTAB GTAB



GCR = GTAB.GCR;
GCR(GCR > .5) = [];

sum(GCR < .001)

%GCR(GCR < .001) = [];


close all; histogram(GCR,0:.001:.1)
[N,E,B] = histcounts(GCR,0:.001:.1);
N = N';
E = E';






VCR = VTAB.VCR;
VCR(VCR > .5) = [];
sum(VCR < .001)
close all; histogram(VCR,0:.0001:.01)
[N,E,B] = histcounts(VCR,0:.0001:.01);
N = N';
E = E';






%% EOF