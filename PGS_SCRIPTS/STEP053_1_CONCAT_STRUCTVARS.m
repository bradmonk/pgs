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
%% LOAD DATASETS
%==========================================================================
clc; clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP052_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;       % GNOMAD DATASET
GUO     = PGS.CLINVAR;      % GUO & GREGG (2019) CLINVAR DATA
MKO     = PGS.MOUSE;        % MOUSE KO DATASET


DAW = load([P.mat P.f 'DAWES_S2_S3.mat']);

SVS = load([P.mat P.f 'STRUCTVARS.mat']);





%==========================================================================
%% RENAME VARIABLES
%==========================================================================
clc; clearvars -except P SVS GNOMAD GUO MKO DAW


STRVAR = SVS.STRUCTVARS;


GNOMAD.Properties.VariableNames{'geneName'} = 'GENE';

DAW.DAWES2_S1.Properties.VariableNames{'gene'} = 'DGENE';

STRVAR.Properties.VariableNames{'x_CHROM'} = 'CHR';





%==========================================================================
%% (STEP 1 of N) CONCATINATE GNOMAD LOF & STRUCT VAR TABLES
%==========================================================================
% 
% Here we start the process of LOF & STRUCT table concatination. In this
% step a CHRPOS variable is made for each table, and the INFO column
% is dropped from the STRVAR table.
% 
%       WITH THE "INFO" COLUMN, THE STRVAR TABLE IS 1.60 GB
%    WITHOUT THE "INFO" COLUMN, THE STRVAR TABLE IS 0.08 GB
% 
% We will also convert all char cells to strings.
% 
%-------------------------
clc; clearvars -except P SVS STRVAR GNOMAD GUO MKO DAW



peek(GNOMAD)
peek(STRVAR)

% TABLE VARIABLE NAMES
%-------------------------
%{
GNOMAD.CHR
GNOMAD.POS
GNOMAD.POS_END
GNOMAD.GENE
GNOMAD.AN
GNOMAD.RC
GNOMAD.AC
GNOMAD.AF
GNOMAD.NHOMALT
GNOMAD.AN_raw
GNOMAD.AN_male
GNOMAD.AN_female
GNOMAD.geneID
GNOMAD.rs
GNOMAD.refBase
GNOMAD.altBase
GNOMAD.altAllele
GNOMAD.N_alleleTypes
GNOMAD.AC_RAW
GNOMAD.AF_RAW
GNOMAD.NHOMALT_RAW
GNOMAD.AC_MALE
GNOMAD.AC_FEMALE
GNOMAD.AF_MALE
GNOMAD.AF_FEMALE
GNOMAD.NHOMALT_MALE
GNOMAD.NHOMALT_FEMALE
GNOMAD.varType
GNOMAD.alleleType
GNOMAD.consequence
GNOMAD.impact
GNOMAD.biotype
GNOMAD.quality
GNOMAD.filter
GNOMAD.rf_tp
GNOMAD.FS
GNOMAD.inbredCoef
GNOMAD.MQ
GNOMAD.Z_mapQuality
GNOMAD.callConf
GNOMAD.Z_posBias
GNOMAD.strandBias
GNOMAD.Z_baseQuality
GNOMAD.Z_hardClip
GNOMAD.depth
GNOMAD.oddsTrueVar
GNOMAD.VQSR_culprit
GNOMAD.inSegDup
GNOMAD.inLowComplex
GNOMAD.inDecoy
GNOMAD.nonpar
GNOMAD.rf_positive_label
GNOMAD.rf_negative_label
GNOMAD.rf_label
GNOMAD.rf_train
GNOMAD.was_mixed
GNOMAD.has_star
GNOMAD.pabMax
GNOMAD.popMax
GNOMAD.popMax_AC
GNOMAD.popMax_AN
GNOMAD.popMax_AF
GNOMAD.popMax_NHOMALT
GNOMAD.featureType
GNOMAD.feature
GNOMAD.exon
GNOMAD.intron
GNOMAD.proteinPos
GNOMAD.aminoAcids
GNOMAD.codons
GNOMAD.existingVar
GNOMAD.alleleNum
GNOMAD.strand
GNOMAD.flag
GNOMAD.varClass
GNOMAD.minimised
GNOMAD.symbolSource
GNOMAD.HGNC_ID
GNOMAD.canonical
GNOMAD.CCDS
GNOMAD.ENSP
GNOMAD.SWISSPROT
GNOMAD.TREMBL
GNOMAD.UNIPARC
GNOMAD.genePheno
GNOMAD.HGVS_offset
GNOMAD.CLIN_SIG
GNOMAD.somatic
GNOMAD.pheno
GNOMAD.lof
GNOMAD.lof_flags
GNOMAD.VCR
GNOMAD.GCR
GNOMAD.GENEi
GNOMAD.N
GNOMAD.P
GNOMAD.R
GNOMAD.A
GNOMAD.AA
GNOMAD.RR
GNOMAD.RA
GNOMAD.fR
GNOMAD.fA
GNOMAD.fRR
GNOMAD.fAA
GNOMAD.fRA
GNOMAD.HW_fRR
GNOMAD.HW_fAA
GNOMAD.HW_fRA


STRVAR.CHR
STRVAR.POS
STRVAR.ID
STRVAR.REF
STRVAR.ALT
STRVAR.QUAL
STRVAR.FILTER
STRVAR.INFO
STRVAR.AN
STRVAR.AC
STRVAR.AF
STRVAR.N_BI_GENOS
STRVAR.N_HOMREF
STRVAR.N_HET
STRVAR.N_HOMALT
STRVAR.FREQ_HOMREF
STRVAR.FREQ_HET
STRVAR.FREQ_HOMALT
%}
%-------------------------




i = isnan(STRVAR.CHR); sum(i)
j = isnan(STRVAR.POS); sum(j)

STRVAR.CHR(i) = 25;
STRVAR.POS(j) = 0;


GNOMAD.CHRPOS = makeCHRPOS(GNOMAD.CHR, GNOMAD.POS);

STRVAR.CHRPOS = makeCHRPOS(STRVAR.CHR, STRVAR.POS);


GNOMAD = movevars(GNOMAD,'CHRPOS','Before','CHR');

STRVAR = movevars(STRVAR,'CHRPOS','Before','CHR');


STRVAR = convertvars(STRVAR,@iscell,'string');



SVS.STRUCTVARS = STRVAR;

STRVAR.INFO = [];




unique(STRVAR.ALT)
unique(STRVAR.CHR)



%==========================================================================
%% (STEP 2 of N) CONCATINATE GNOMAD LOF & STRUCT VAR TABLES
%==========================================================================
% 
% The "GNOMAD" LoF table has many more variables (128989 x 110) than the
% "STRVAR" structural variant table (217308 x 18), but less sites.
% 
% Here we start to match variable columns in GNOMAD and STRVAR.
% 
%{
% 
% At single locus with two alleles R and A with frequencies fR and fA,
% respectively, the expected genotype frequencies under random mating are:
% 
% fRR = fR^2
% fAA = fA^2
% fRA = 2*fR*fA
%
% N.B. HARDY-WEINBERG DOES NOT WORK FOR SEX CHROMOSOMES


GNO = GNOMAD;

GNO.RC = GNO.AN - GNO.AC;


% GIVEN...
N   = GNO.AN;
A   = GNO.AC;
R   = GNO.RC;
AA  = GNO.NHOMALT;



% COMPUTED
RA          = A - (AA .* 2);
RR          = (R - RA) ./ 2;
RR_AA_RA    = RR + AA + RA;
R_A         = R + A;
PPL         = N ./ 2;
fR          = R ./ N;
fA          = A ./ N;
fRR         = RR ./ PPL;
fAA         = AA ./ PPL;
fRA         = RA ./ PPL;

% HARDY-WEINBERG ESTIMATES
HW_fRR      = fR.^2;
HW_fAA      = fA.^2;
HW_fRA      = 2 .* fR .* fA;


FREQ = table(N, R, A, AA, RR, RA, R_A , RR_AA_RA, PPL,...
             fR, fA, fRR, fAA, fRA, HW_fRR, HW_fAA, HW_fRA);


GNO.N   = N;
GNO.P   = PPL;
GNO.R   = R;
GNO.A   = A;
GNO.AA  = AA;
GNO.RR  = RR;
GNO.RA  = RA;
GNO.fR  = fR;
GNO.fA  = fA;
GNO.fRR = fRR;
GNO.fAA = fAA;
GNO.fRA = fRA;
GNO.HW_fRR = HW_fRR;
GNO.HW_fAA = HW_fAA;
GNO.HW_fRA = HW_fRA;
%}
%-------------------------
clc; clearvars -except P SVS STRVAR GNOMAD GUO MKO DAW


GNOMAD.isLOF = ones(size(GNOMAD,1),1);
GNOMAD.isSV  = zeros(size(GNOMAD,1),1);


SVARS = GNOMAD;
SVARS(1:end,:) = [];
SVARS.isSV(size(STRVAR,1)) = 1;



SVARS.CHRPOS        = STRVAR.CHRPOS;
SVARS.CHR           = STRVAR.CHR;
SVARS.POS           = STRVAR.POS;
SVARS.refBase       = STRVAR.REF;
SVARS.altBase       = STRVAR.ALT;
SVARS.altAllele     = STRVAR.ALT;
SVARS.varType       = STRVAR.ALT;
SVARS.alleleType    = STRVAR.ALT;
SVARS.consequence   = STRVAR.ID;
SVARS.filter        = STRVAR.FILTER;
SVARS.isLOF         = zeros(size(STRVAR,1),1);
SVARS.isSV          = ones(size(STRVAR,1),1);
SVARS.quality       = STRVAR.QUAL;
SVARS.AN            = STRVAR.AN;
SVARS.AC            = STRVAR.AC;
SVARS.RC            = SVARS.AN - SVARS.AC;
SVARS.AF            = STRVAR.AF;
SVARS.NHOMALT       = STRVAR.N_HOMALT;
SVARS.N             = STRVAR.AN;
SVARS.P             = STRVAR.AN ./ 2;
SVARS.R             = SVARS.RC;
SVARS.A             = SVARS.AC;
SVARS.AA            = STRVAR.N_HOMALT;
SVARS.RA            = SVARS.A - (SVARS.AA .* 2);
SVARS.RR            = (SVARS.R - SVARS.RA) ./ 2;
SVARS.fR            = SVARS.R ./ SVARS.N;
SVARS.fA            = SVARS.A ./ SVARS.N;
SVARS.fRR           = SVARS.RR ./ SVARS.P;
SVARS.fAA           = SVARS.AA ./ SVARS.P;
SVARS.fRA           = SVARS.RA ./ SVARS.P;
SVARS.HW_fRR        = SVARS.fR.^2;
SVARS.HW_fAA        = SVARS.fA.^2;
SVARS.HW_fRA        = 2 .* SVARS.fR .* SVARS.fA;



% VARIABLES IN STRVAR TABLE BUT NOT IN GNOMAD TABLE
% SVARS.N_BI_GENOS    = STRVAR.N_BI_GENOS;    % REDUNDANT WITH SVARS.P
% SVARS.N_HOMREF      = STRVAR.N_HOMREF;      % REDUNDANT WITH SVARS.RR
% SVARS.N_HET         = STRVAR.N_HET;         % REDUNDANT WITH SVARS.RA
% SVARS.FREQ_HOMREF   = STRVAR.FREQ_HOMREF;   % REDUNDANT WITH SVARS.fRR
% SVARS.FREQ_HET      = STRVAR.FREQ_HET;      % REDUNDANT WITH SVARS.fRA
% SVARS.FREQ_HOMALT   = STRVAR.FREQ_HOMALT;   % SHOULD BE REDUNDANT WITH SVARS.fAA (but isnt ??)


GNOM = GNOMAD;

GNOMAD = [GNOMAD; SVARS];





%==========================================================================
%% (STEP 3 of N) SORT GNOMAD TABLE AND FILL MISSING GENES WITH PREV
%==========================================================================
clc; clearvars -except P SVS STRVAR SVARS GNOMAD GNOM GUO MKO DAW


GNO = GNOMAD;

GNO = sortrows(GNO,'CHRPOS');

GNO.GENE(1) = "OR4F5";

GNO.GENE = fillmissing(GNO.GENE,'previous');

GNOMAD = GNO;





%==========================================================================
%% (STEP 4 of N) MATCH THE QUALITY SCORE COLUMNS IN GNOMAD & STRVAR
%==========================================================================
% 
% The "GNOMAD" LoF table has many more variables (128989 x 110) than the
% "STRVAR" structural variant table (217308 x 18), but less sites.
% 
% Here we start to match variable columns in GNOMAD and STRVAR.
% 
%-------------------------
%{
clc; clearvars -except P SVS STRVAR SVARS GNOMAD GNOM GUO MKO DAW


GNOMAD_logQUAL = log(GNOMAD.quality);
STRVAR_logQUAL = log(STRVAR.QUAL);

close all;
histogram(GNOMAD_logQUAL)
histogram(STRVAR_logQUAL)


GNOMAD_logQUAL = uint64(rescale(GNOMAD_logQUAL,0,100));
STRVAR_logQUAL = uint64(rescale(STRVAR_logQUAL,35,100));


figure
subplot(1,2,1); histogram(GNOMAD_logQUAL)
subplot(1,2,2); histogram(STRVAR_logQUAL)


figure
histogram(GNOMAD_logQUAL)
hold on;
histogram(STRVAR_logQUAL)
%}








%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P SVS STRVAR SVARS GNOMAD GNOM GUO MKO DAW


save([P.mat P.f 'PGS_STEP053_1_CONCAT_STRUCTVARS.mat'],'GNOMAD','GNOM','MKO','GUO','DAW');


% writetable(GNOMAD,[P.csv P.f 'PGS_STEP054_GNOMAD.csv'])
% writetable(DMAD,[P.csv P.f 'PGS_STEP054_DMAD.csv'])
% writetable(DAWMAD,[P.csv P.f 'PGS_STEP055_DAWMAD.csv'])
