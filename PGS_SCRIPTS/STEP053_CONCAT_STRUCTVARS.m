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



% SEE BLOG POST ABOUT GNOMAD STRUCTURAL VARIANTS:
% https://macarthurlab.org/2019/03/20/structural-variants-in-gnomad/
% https://jmonlong.github.io/Hippocamplus/2019/03/31/gnomad-sv-first-look/


%==========================================================================
%% LOAD DATASETS
%==========================================================================
clc; clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP052_OUTPUT.mat']);






%==========================================================================
%% CONCATINATE GNOMAD & SVLOF
%==========================================================================
clc; clearvars -except P PGS


GNOMAD = PGS.GNOMAD;
SVLOF  = PGS.SVLOF;



GNO = [GNOMAD; SVLOF];


i = isnan(GNO.CHR); sum(i)
j = isnan(GNO.POS); sum(j)
GNO.CHR(i) = 26;
GNO.POS(j) = 0;
GNO.CHRPOS = makeCHRPOS(GNO.CHR, GNO.POS);


GNO = sortrows(GNO,'CHRPOS');


GNO = GENEijk(GNO);






%==========================================================================
%% ORGANIZE TABLE VARS
%==========================================================================
clc; clearvars -except P PGS GNO



TBL = GNO;


clc; disp(char(string(TBL.Properties.VariableNames')));





TBL.faf95       = [];
TBL.faf95_afr   = [];
TBL.faf95_amr   = [];
TBL.faf95_eas   = [];
TBL.faf95_nfe   = [];
TBL.faf95_sas   = [];
TBL.faf99       = [];
TBL.faf99_afr   = [];
TBL.faf99_amr   = [];
TBL.faf99_eas   = [];
TBL.faf99_nfe   = [];
TBL.faf99_sas   = [];




TBL = movevars(TBL,{'n_alt_alleles'},'After','ALLELE_NUM');

TBL = movevars(TBL,{'nhomalt'},'After','AC');

TBL = movevars(TBL,  {...
'nhomalt_raw',...
'nhomalt_popmax',...
'nhomalt_female',...
'nhomalt_male',...
'nhomalt_afr',...
'nhomalt_amr',...
'nhomalt_asj',...
'nhomalt_eas',...
'nhomalt_eas_jpn',...
'nhomalt_eas_kor',...
'nhomalt_eas_oea',...
'nhomalt_fin',...
'nhomalt_nfe',...
'nhomalt_nfe_bgr',...
'nhomalt_nfe_est',...
'nhomalt_nfe_nwe',...
'nhomalt_nfe_onf',...
'nhomalt_nfe_seu',...
'nhomalt_nfe_swe',...
'nhomalt_oth',...
'nhomalt_sas',...
},'After','AC_sas');






TBL = movevars(TBL,  {...
'N'
'P'
'R'
'A'
'AA'
'RR'
'RA'
'fR'
'fA'
'fRR'
'fAA'
'fRA'
'HW_fRR'
'HW_fAA'
'HW_fRA'
},'After','nhomalt');






clc; disp(char(string(TBL.Properties.VariableNames')));



GNO = TBL;



%==========================================================================
%% COMPUTE VCR & GCR
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
% 
% 
% 
% 
% 
% GCR(g) = 1 - prod(1 - VCR(1:k));
% 
% 
% The equation above computes the gene carrier rate (GCR) for a gene 'g',
% where VCR is the variant carrier rates. The product is computed for all 
% variants 1:k of a given gene.
% 
%--------------------------------------------------------------------------
clc; clearvars -except P PGS GNO


TBL = GNO;


TBL.VCR           = vcr(TBL.AN,         TBL.AC,         TBL.nhomalt);
TBL.VCR_raw       = vcr(TBL.AN_raw,     TBL.AC_raw,     TBL.nhomalt_raw);
TBL.VCR_popmax    = vcr(TBL.AN_popmax,  TBL.AC_popmax,  TBL.nhomalt_popmax);
TBL.VCR_afr       = vcr(TBL.AN_afr,     TBL.AC_afr,     TBL.nhomalt_afr);
TBL.VCR_amr       = vcr(TBL.AN_amr,     TBL.AC_amr,     TBL.nhomalt_amr);
TBL.VCR_asj       = vcr(TBL.AN_asj,     TBL.AC_asj,     TBL.nhomalt_asj);
TBL.VCR_eas       = vcr(TBL.AN_eas,     TBL.AC_eas,     TBL.nhomalt_eas);
TBL.VCR_fin       = vcr(TBL.AN_fin,     TBL.AC_fin,     TBL.nhomalt_fin);
TBL.VCR_nfe       = vcr(TBL.AN_nfe,     TBL.AC_nfe,     TBL.nhomalt_nfe);
TBL.VCR_oth       = vcr(TBL.AN_oth,     TBL.AC_oth,     TBL.nhomalt_oth);
TBL.VCR_sas       = vcr(TBL.AN_sas,     TBL.AC_sas,     TBL.nhomalt_sas);
TBL.VCR_eas_jpn   = vcr(TBL.AN_eas_jpn, TBL.AC_eas_jpn, TBL.nhomalt_eas_jpn);
TBL.VCR_eas_kor   = vcr(TBL.AN_eas_kor, TBL.AC_eas_kor, TBL.nhomalt_eas_kor);
TBL.VCR_eas_oea   = vcr(TBL.AN_eas_oea, TBL.AC_eas_oea, TBL.nhomalt_eas_oea);
TBL.VCR_nfe_bgr   = vcr(TBL.AN_nfe_bgr, TBL.AC_nfe_bgr, TBL.nhomalt_nfe_bgr);
TBL.VCR_nfe_est   = vcr(TBL.AN_nfe_est, TBL.AC_nfe_est, TBL.nhomalt_nfe_est);
TBL.VCR_nfe_nwe   = vcr(TBL.AN_nfe_nwe, TBL.AC_nfe_nwe, TBL.nhomalt_nfe_nwe);
TBL.VCR_nfe_onf   = vcr(TBL.AN_nfe_onf, TBL.AC_nfe_onf, TBL.nhomalt_nfe_onf);
TBL.VCR_nfe_seu   = vcr(TBL.AN_nfe_seu, TBL.AC_nfe_seu, TBL.nhomalt_nfe_seu);
TBL.VCR_nfe_swe   = vcr(TBL.AN_nfe_swe, TBL.AC_nfe_swe, TBL.nhomalt_nfe_swe);




TBL.GCR           = gcr(TBL.VCR,         TBL.GENEi);
TBL.GCR_raw       = gcr(TBL.VCR_raw,     TBL.GENEi);
TBL.GCR_popmax    = gcr(TBL.VCR_popmax,  TBL.GENEi);
TBL.GCR_afr       = gcr(TBL.VCR_afr,     TBL.GENEi);
TBL.GCR_amr       = gcr(TBL.VCR_amr,     TBL.GENEi);
TBL.GCR_asj       = gcr(TBL.VCR_asj,     TBL.GENEi);
TBL.GCR_eas       = gcr(TBL.VCR_eas,     TBL.GENEi);
TBL.GCR_fin       = gcr(TBL.VCR_fin,     TBL.GENEi);
TBL.GCR_nfe       = gcr(TBL.VCR_nfe,     TBL.GENEi);
TBL.GCR_oth       = gcr(TBL.VCR_oth,     TBL.GENEi);
TBL.GCR_sas       = gcr(TBL.VCR_sas,     TBL.GENEi);
TBL.GCR_eas_jpn   = gcr(TBL.VCR_eas_jpn, TBL.GENEi);
TBL.GCR_eas_kor   = gcr(TBL.VCR_eas_kor, TBL.GENEi);
TBL.GCR_eas_oea   = gcr(TBL.VCR_eas_oea, TBL.GENEi);
TBL.GCR_nfe_bgr   = gcr(TBL.VCR_nfe_bgr, TBL.GENEi);
TBL.GCR_nfe_est   = gcr(TBL.VCR_nfe_est, TBL.GENEi);
TBL.GCR_nfe_nwe   = gcr(TBL.VCR_nfe_nwe, TBL.GENEi);
TBL.GCR_nfe_onf   = gcr(TBL.VCR_nfe_onf, TBL.GENEi);
TBL.GCR_nfe_seu   = gcr(TBL.VCR_nfe_seu, TBL.GENEi);
TBL.GCR_nfe_swe   = gcr(TBL.VCR_nfe_swe, TBL.GENEi);



GNO = TBL;




%==========================================================================
%% (STEP 3 of N) SORT GNOMAD TABLE AND FILL MISSING GENES WITH PREV
%==========================================================================
% clc; clearvars -except P PGS GNO
% 
% 
% TBL = GNO;
% 
% TBL = sortrows(TBL,'CHRPOS');
% 
% TBL.GENE = fillmissing(TBL.GENE,'previous');
% 
% GNO = TBL;
%}
%-------------------------











%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNO



GNOMAD          = GNO;
DAWES           = PGS.DAWES;
DATAS           = PGS.DATAS;
DATAS.SVFULL    = PGS.SVFULL;
DATAS.SVLOF     = PGS.SVLOF;


save([P.mat P.f 'PGS_STEP053_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');



% DATAS = PGS.DATAS;
% DAWES = PGS.DAWES;
% save([P.mat P.f 'PGS_STEP052_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');



