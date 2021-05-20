%==========================================================================
%% SCRIPT_NAME
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
%------------------------------------------------





%==========================================================================
%% LOAD DATASET
%==========================================================================
clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP030_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;







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
clc; clearvars -except P PGS GNOMAD


peek(GNOMAD)


GNOMAD.VCR = (GNOMAD.AC - GNOMAD.nhomalt) ./ (.5 .* GNOMAD.AN);

GNOMAD.VCR_raw          = (GNOMAD.AC_raw - GNOMAD.nhomalt_raw)          ./ (.5 .* GNOMAD.AN_raw);
GNOMAD.VCR_popmax       = (GNOMAD.AC_popmax - GNOMAD.nhomalt_popmax)    ./ (.5 .* GNOMAD.AN_popmax);

GNOMAD.VCR_afr          = (GNOMAD.AC_afr - GNOMAD.nhomalt_afr)          ./ (.5 .* GNOMAD.AN_afr);
GNOMAD.VCR_amr          = (GNOMAD.AC_amr - GNOMAD.nhomalt_amr)          ./ (.5 .* GNOMAD.AN_amr);
GNOMAD.VCR_asj          = (GNOMAD.AC_asj - GNOMAD.nhomalt_asj)          ./ (.5 .* GNOMAD.AN_asj);
GNOMAD.VCR_eas          = (GNOMAD.AC_eas - GNOMAD.nhomalt_eas)          ./ (.5 .* GNOMAD.AN_eas);
GNOMAD.VCR_eas_jpn      = (GNOMAD.AC_eas_jpn - GNOMAD.nhomalt_eas_jpn)  ./ (.5 .* GNOMAD.AN_eas_jpn);
GNOMAD.VCR_eas_kor      = (GNOMAD.AC_eas_kor - GNOMAD.nhomalt_eas_kor)  ./ (.5 .* GNOMAD.AN_eas_kor);
GNOMAD.VCR_eas_oea      = (GNOMAD.AC_eas_oea - GNOMAD.nhomalt_eas_oea)  ./ (.5 .* GNOMAD.AN_eas_oea);
GNOMAD.VCR_fin          = (GNOMAD.AC_fin - GNOMAD.nhomalt_fin)          ./ (.5 .* GNOMAD.AN_fin);
GNOMAD.VCR_nfe          = (GNOMAD.AC_nfe - GNOMAD.nhomalt_nfe)          ./ (.5 .* GNOMAD.AN_nfe);
GNOMAD.VCR_nfe_bgr      = (GNOMAD.AC_nfe_bgr - GNOMAD.nhomalt_nfe_bgr)  ./ (.5 .* GNOMAD.AN_nfe_bgr);
GNOMAD.VCR_nfe_est      = (GNOMAD.AC_nfe_est - GNOMAD.nhomalt_nfe_est)  ./ (.5 .* GNOMAD.AN_nfe_est);
GNOMAD.VCR_nfe_nwe      = (GNOMAD.AC_nfe_nwe - GNOMAD.nhomalt_nfe_nwe)  ./ (.5 .* GNOMAD.AN_nfe_nwe);
GNOMAD.VCR_nfe_onf      = (GNOMAD.AC_nfe_onf - GNOMAD.nhomalt_nfe_onf)  ./ (.5 .* GNOMAD.AN_nfe_onf);
GNOMAD.VCR_nfe_seu      = (GNOMAD.AC_nfe_seu - GNOMAD.nhomalt_nfe_seu)  ./ (.5 .* GNOMAD.AN_nfe_seu);
GNOMAD.VCR_nfe_swe      = (GNOMAD.AC_nfe_swe - GNOMAD.nhomalt_nfe_swe)  ./ (.5 .* GNOMAD.AN_nfe_swe);
GNOMAD.VCR_oth          = (GNOMAD.AC_oth - GNOMAD.nhomalt_oth)          ./ (.5 .* GNOMAD.AN_oth);
GNOMAD.VCR_sas          = (GNOMAD.AC_sas - GNOMAD.nhomalt_sas)          ./ (.5 .* GNOMAD.AN_sas);









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
clc; clearvars -except P PGS GNOMAD







% ADD GCR COLUMN TO TABLE
GNOMAD.GCR = nan(height(GNOMAD),1);

GNOMAD.GCR_raw      = nan(height(GNOMAD),1);
GNOMAD.GCR_popmax   = nan(height(GNOMAD),1);
GNOMAD.GCR_afr      = nan(height(GNOMAD),1);
GNOMAD.GCR_amr      = nan(height(GNOMAD),1);
GNOMAD.GCR_asj      = nan(height(GNOMAD),1);
GNOMAD.GCR_eas      = nan(height(GNOMAD),1);
GNOMAD.GCR_eas_jpn  = nan(height(GNOMAD),1);
GNOMAD.GCR_eas_kor  = nan(height(GNOMAD),1);
GNOMAD.GCR_eas_oea  = nan(height(GNOMAD),1);
GNOMAD.GCR_fin      = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe      = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe_bgr  = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe_est  = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe_nwe  = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe_onf  = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe_seu  = nan(height(GNOMAD),1);
GNOMAD.GCR_nfe_swe  = nan(height(GNOMAD),1);
GNOMAD.GCR_oth      = nan(height(GNOMAD),1);
GNOMAD.GCR_sas      = nan(height(GNOMAD),1);


GNOMAD = GENEijk(GNOMAD);


% COMPUTE GCR FOR EACH GENE
for g = 1:numel(unique(GNOMAD.GENEi))

    GNOMAD.GCR(GNOMAD.GENEi==g) = 1 - prod(1 - GNOMAD.VCR(GNOMAD.GENEi==g));

    GNOMAD.GCR_raw(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_raw(GNOMAD.GENEi==g));
    GNOMAD.GCR_popmax(GNOMAD.GENEi==g)       = 1 - prod(1 - GNOMAD.VCR_popmax(GNOMAD.GENEi==g));
    GNOMAD.GCR_afr(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_afr(GNOMAD.GENEi==g));
    GNOMAD.GCR_amr(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_amr(GNOMAD.GENEi==g));
    GNOMAD.GCR_asj(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_asj(GNOMAD.GENEi==g));
    GNOMAD.GCR_eas(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_eas(GNOMAD.GENEi==g));
    GNOMAD.GCR_eas_jpn(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_eas_jpn(GNOMAD.GENEi==g));
    GNOMAD.GCR_eas_kor(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_eas_kor(GNOMAD.GENEi==g));
    GNOMAD.GCR_eas_oea(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_eas_oea(GNOMAD.GENEi==g));
    GNOMAD.GCR_fin(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_fin(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_nfe(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe_bgr(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_nfe_bgr(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe_est(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_nfe_est(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe_nwe(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_nfe_nwe(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe_onf(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_nfe_onf(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe_seu(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_nfe_seu(GNOMAD.GENEi==g));
    GNOMAD.GCR_nfe_swe(GNOMAD.GENEi==g)      = 1 - prod(1 - GNOMAD.VCR_nfe_swe(GNOMAD.GENEi==g));
    GNOMAD.GCR_oth(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_oth(GNOMAD.GENEi==g));
    GNOMAD.GCR_sas(GNOMAD.GENEi==g)          = 1 - prod(1 - GNOMAD.VCR_sas(GNOMAD.GENEi==g));


end








%==========================================================================
%% COMPUTE CUMULATIVE CARRIER RATE (CCR)
%==========================================================================
% 
% CCR = 1 - prod(1 - GCR(L) )
% 
% 
% The equation above cumulative carrier (CCR) for a list of genes 'L',
% where GCR is an array of gene carrier rates. The product is computed 
% for the 1-GCR values of all genes in list L.
% 
%--------------------------------------------------------------------------
% % EXAMPLE...
% 
% 
% n = height(GEN);
% 
% GCR = GEN.GCR;
% 
% 
% cc = 1000;
% gg = 30;
% 
% CCR = nan(gg,cc);
% 
% for c = 1:cc
% for g = 1:gg
% 
%     r = randperm(n,g);
% 
%     CCR(g,c) = 1 - prod(1 - GCR(r));
% 
% end
% end
% 
% 
%--------------------------------------------------------------------------











%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DATAS = PGS.DATAS;


save([P.mat P.f 'PGS_STEP040_OUTPUT.mat'],'GNOMAD','DATAS');



% writetable(GNOMAD,[P.csv P.f 'PGS_STEP040_OUTPUT.csv'])














