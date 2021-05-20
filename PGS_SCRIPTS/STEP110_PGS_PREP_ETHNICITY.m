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
P.fun    = [P.home filesep 'PGS_FUNS'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------







%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P



PGS = load([P.mat P.f 'PGS_STEP100_OUTPUT.mat']);


GNOMAD  = PGS.GNOMAD; peek(GNOMAD)

GNOMAD = GENEijk(GNOMAD);





peek(GNOMAD); gene_report(GNOMAD)
%==========================================================================
%% PREPROCESS GNOMAD
%==========================================================================
clc; clearvars -except P PGS GNOMAD




% % REMOVE VARIANTS THAT HAVE OVER X% MAF
% %------------------------------------------------------------------
% GNOMAD(GNOMAD.AF > .50 , :) = [];
% GNOMAD = GENEijk(GNOMAD);
% 
% 
% 
% 
% 
% % REMOVE VARIANTS THAT HAVE MORE THAN N% HOMOZYGOUS MINOR ALLELES
% %------------------------------------------------------------------
% pct_nhomalt = (GNOMAD.nhomalt ./ GNOMAD.AC);
% i = pct_nhomalt > .05;
% GNOMAD(i , :) = [];
% GNOMAD = GENEijk(GNOMAD);










%==========================================================================
%% COMPUTE ETHNICITY VCR
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
clc; clearvars -except P PGS GNOMAD



GNO = GNOMAD;



GNO.VCR              = (GNO.AC - GNO.nhomalt) ./ (.5 .* GNO.AN);

GNO.VCR_raw          = (GNO.AC_raw - GNO.nhomalt_raw)          ./ (.5 .* GNO.AN_raw);
GNO.VCR_popmax       = (GNO.AC_popmax - GNO.nhomalt_popmax)    ./ (.5 .* GNO.AN_popmax);
GNO.VCR_male         = (GNO.AC_male - GNO.nhomalt_male)        ./ (.5 .* GNO.AN_male);
GNO.VCR_female       = (GNO.AC_female - GNO.nhomalt_female)    ./ (.5 .* GNO.AN_female);

GNO.VCR_afr          = (GNO.AC_afr - GNO.nhomalt_afr)          ./ (.5 .* GNO.AN_afr);
GNO.VCR_amr          = (GNO.AC_amr - GNO.nhomalt_amr)          ./ (.5 .* GNO.AN_amr);
GNO.VCR_asj          = (GNO.AC_asj - GNO.nhomalt_asj)          ./ (.5 .* GNO.AN_asj);
GNO.VCR_eas          = (GNO.AC_eas - GNO.nhomalt_eas)          ./ (.5 .* GNO.AN_eas);
GNO.VCR_eas_jpn      = (GNO.AC_eas_jpn - GNO.nhomalt_eas_jpn)  ./ (.5 .* GNO.AN_eas_jpn);
GNO.VCR_eas_kor      = (GNO.AC_eas_kor - GNO.nhomalt_eas_kor)  ./ (.5 .* GNO.AN_eas_kor);
GNO.VCR_eas_oea      = (GNO.AC_eas_oea - GNO.nhomalt_eas_oea)  ./ (.5 .* GNO.AN_eas_oea);
GNO.VCR_fin          = (GNO.AC_fin - GNO.nhomalt_fin)          ./ (.5 .* GNO.AN_fin);
GNO.VCR_nfe          = (GNO.AC_nfe - GNO.nhomalt_nfe)          ./ (.5 .* GNO.AN_nfe);
GNO.VCR_nfe_bgr      = (GNO.AC_nfe_bgr - GNO.nhomalt_nfe_bgr)  ./ (.5 .* GNO.AN_nfe_bgr);
GNO.VCR_nfe_est      = (GNO.AC_nfe_est - GNO.nhomalt_nfe_est)  ./ (.5 .* GNO.AN_nfe_est);
GNO.VCR_nfe_nwe      = (GNO.AC_nfe_nwe - GNO.nhomalt_nfe_nwe)  ./ (.5 .* GNO.AN_nfe_nwe);
GNO.VCR_nfe_onf      = (GNO.AC_nfe_onf - GNO.nhomalt_nfe_onf)  ./ (.5 .* GNO.AN_nfe_onf);
GNO.VCR_nfe_seu      = (GNO.AC_nfe_seu - GNO.nhomalt_nfe_seu)  ./ (.5 .* GNO.AN_nfe_seu);
GNO.VCR_nfe_swe      = (GNO.AC_nfe_swe - GNO.nhomalt_nfe_swe)  ./ (.5 .* GNO.AN_nfe_swe);
GNO.VCR_oth          = (GNO.AC_oth - GNO.nhomalt_oth)          ./ (.5 .* GNO.AN_oth);
GNO.VCR_sas          = (GNO.AC_sas - GNO.nhomalt_sas)          ./ (.5 .* GNO.AN_sas);




% VARS = {'VCR','VCR_raw','VCR_popmax','VCR_male','VCR_female','VCR_afr',...
%         'VCR_amr','VCR_asj','VCR_eas','VCR_eas_jpn','VCR_eas_kor',...
%         'VCR_eas_oea','VCR_fin','VCR_nfe','VCR_nfe_bgr','VCR_nfe_est',...
%         'VCR_nfe_nwe','VCR_nfe_onf','VCR_nfe_seu','VCR_nfe_swe',...
%         'VCR_oth','VCR_sas'};
% 
% GNO = movevars(GNO,VARS,'After','AF');


GNOMAD = GNO;



%==========================================================================
%% COMPUTE ETHNICITY GCR
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
clc; clearvars -except P PGS GNOMAD



GNOMAD = GENEijk(GNOMAD);

GNO = GNOMAD;



% ADD GCR COLUMN TO TABLE
GNO.GCR          = nan(height(GNO),1);

GNO.GCR_raw      = nan(height(GNO),1);
GNO.GCR_popmax   = nan(height(GNO),1);
GNO.GCR_male     = nan(height(GNO),1);
GNO.GCR_female   = nan(height(GNO),1);

GNO.GCR_afr      = nan(height(GNO),1);
GNO.GCR_amr      = nan(height(GNO),1);
GNO.GCR_asj      = nan(height(GNO),1);
GNO.GCR_eas      = nan(height(GNO),1);
GNO.GCR_eas_jpn  = nan(height(GNO),1);
GNO.GCR_eas_kor  = nan(height(GNO),1);
GNO.GCR_eas_oea  = nan(height(GNO),1);
GNO.GCR_fin      = nan(height(GNO),1);
GNO.GCR_nfe      = nan(height(GNO),1);
GNO.GCR_nfe_bgr  = nan(height(GNO),1);
GNO.GCR_nfe_est  = nan(height(GNO),1);
GNO.GCR_nfe_nwe  = nan(height(GNO),1);
GNO.GCR_nfe_onf  = nan(height(GNO),1);
GNO.GCR_nfe_seu  = nan(height(GNO),1);
GNO.GCR_nfe_swe  = nan(height(GNO),1);
GNO.GCR_oth      = nan(height(GNO),1);
GNO.GCR_sas      = nan(height(GNO),1);




% COMPUTE GCR FOR EACH GENE
for g = 1:numel(unique(GNO.GENEi))

    GNO.GCR(GNO.GENEi==g)           = 1 - prod(1 - GNO.VCR(GNO.GENEi==g));

    GNO.GCR_raw(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_raw(GNO.GENEi==g));
    GNO.GCR_popmax(GNO.GENEi==g)    = 1 - prod(1 - GNO.VCR_popmax(GNO.GENEi==g));
    GNO.GCR_male(GNO.GENEi==g)      = 1 - prod(1 - GNO.VCR_male(GNO.GENEi==g));
    GNO.GCR_female(GNO.GENEi==g)    = 1 - prod(1 - GNO.VCR_female(GNO.GENEi==g));

    GNO.GCR_afr(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_afr(GNO.GENEi==g));
    GNO.GCR_amr(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_amr(GNO.GENEi==g));
    GNO.GCR_asj(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_asj(GNO.GENEi==g));
    GNO.GCR_eas(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_eas(GNO.GENEi==g));
    GNO.GCR_eas_jpn(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_eas_jpn(GNO.GENEi==g));
    GNO.GCR_eas_kor(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_eas_kor(GNO.GENEi==g));
    GNO.GCR_eas_oea(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_eas_oea(GNO.GENEi==g));
    GNO.GCR_fin(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_fin(GNO.GENEi==g));
    GNO.GCR_nfe(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_nfe(GNO.GENEi==g));
    GNO.GCR_nfe_bgr(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_nfe_bgr(GNO.GENEi==g));
    GNO.GCR_nfe_est(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_nfe_est(GNO.GENEi==g));
    GNO.GCR_nfe_nwe(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_nfe_nwe(GNO.GENEi==g));
    GNO.GCR_nfe_onf(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_nfe_onf(GNO.GENEi==g));
    GNO.GCR_nfe_seu(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_nfe_seu(GNO.GENEi==g));
    GNO.GCR_nfe_swe(GNO.GENEi==g)   = 1 - prod(1 - GNO.VCR_nfe_swe(GNO.GENEi==g));
    GNO.GCR_oth(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_oth(GNO.GENEi==g));
    GNO.GCR_sas(GNO.GENEi==g)       = 1 - prod(1 - GNO.VCR_sas(GNO.GENEi==g));


end




% VARS = {'GCR','GCR_raw','GCR_popmax','GCR_male','GCR_female','GCR_afr',...
%         'GCR_amr','GCR_asj','GCR_eas','GCR_eas_jpn','GCR_eas_kor',...
%         'GCR_eas_oea','GCR_fin','GCR_nfe','GCR_nfe_bgr','GCR_nfe_est',...
%         'GCR_nfe_nwe','GCR_nfe_onf','GCR_nfe_seu','GCR_nfe_swe',...
%         'GCR_oth','GCR_sas'};
% 
% GNO = movevars(GNO,VARS,'After','VCR_sas');


GNOMAD = GNO;



%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DATAS = PGS.DATAS;

DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP110_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');



