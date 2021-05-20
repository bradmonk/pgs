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


PGS = load([P.mat P.f 'PGS_STEP048_OUTPUT.mat']);


SVS = load([P.mat P.f 'STRUCTVARS.mat']);


% RENAME SVS FIELDS
SVDATA  = SVS.STRUCTVARS;
SVINFO  = SVS.SVARS;
clear SVS
SVS.SVDATA = SVDATA;
SVS.SVINFO = SVINFO;




%==========================================================================
%% CLEAN SVDATA TABLE
%==========================================================================
clc; clearvars -except P PGS SVS



SVDATA = SVS.SVDATA;


SVDATA.Properties.VariableNames{'x_CHROM'} = 'CHR';


SVDATA = convertvars(SVDATA,@iscell,'string');


SVDATA.CHRPOS = makeCHRPOS(SVDATA.CHR,SVDATA.POS);

 
SVDATA = movevars(SVDATA,{'CHRPOS'},'Before','CHR');


SVS.SVDATA = SVDATA;



%==========================================================================
%% LIST THE TYPES OF STRUCTURAL VARIANTS
%==========================================================================
clc; clearvars -except P PGS SVS


SVDATA = SVS.SVDATA;

gcounts( SVDATA.ALT )



% Group           Count 
% ----------      ---------------
% <DEL>            110826 (51.0%)   deletion
% <INS:ME:ALU>      44137 (20.3%)   insertion mobile element Alu
% <DUP>             31029 (14.3%)   duplication
% <INS>             15161 (7.0%)    insertion novel sequence
% <INS:ME:LINE1>     7167 (3.3%)    insertion mobile element LINE1
% <INS:ME:SVA>       4154 (1.9%)    insertion mobile element SINE or VNTR or Alu
% <CPX>              3807 (1.8%)    complex structural variant
% <INV>               573 (0.3%)    inversion
% <INS:ME>            447 (0.2%)    insertion mobile element unspecified
% <CTX>                 7 (0.0%)    translocation



% SV types
% ----------
%   (DEL)   deletion 
%   (INS)   insertion
%   (DUP)   duplication
%   (INV)   inversion 
%   (CTX)   translocation 
%   (CPX)   complex 
%   (ALU)   Alu - aka transposable element
%   (VNTR)  variable number tandem repeat
%   (LINE1) long interspersed nuclear element-1 
%   (SINE)  short interspersed nuclear elements
%   (SVA)   SINE-VNTR-Alu
%   (ME)    mobile element


% ***Alu***
% Alu elements are among the most successful of all mobile elements,
% having over 1 million copies in the human genome (11% of the genome).
% They belong to a class of primate-specific retroelements termed 
% ***SINEs***.
% Short Interspersed Elements (SINEs) are not autonomous, since they 
% acquire trans-acting factors for their amplification from the only 
% active family of autonomous human retroelements:
% ***LINE-1***.


% FUN FACT - Alu etymology
%   Alu elements were originally characterized using the restriction 
%   endonuclease from gram-positive bacteria Arthrobacter luteus (Alu) 


% ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
% ##ALT=<ID=DEL,Description="Deletion">
% ##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">
% ##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">
% ##ALT=<ID=DUP,Description="Duplication">
% ##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
% ##ALT=<ID=INS,Description="Insertion of novel sequence">
% ##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
% ##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
% ##ALT=<ID=INV,Description="Inversion">
% ##ALT=<ID=CNV,Description="Copy number variable region">
% 
% ##ALT=<ID=BND,Description="Translocation">
% ##ALT=<ID=CPX,Description="Complex SV">
% ##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">
% ##ALT=<ID=DEL,Description="Deletion">
% ##ALT=<ID=DUP,Description="Duplication">
% ##ALT=<ID=INS,Description="Insertion">
% ##ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">
% ##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">
% ##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">
% ##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">
% ##ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">
% ##ALT=<ID=INV,Description="Inversion">
% ##FILTER=<ID=MULTIALLELIC,Description="Multiallelic site">



%==========================================================================
%% ONLY USE STRUCTURAL DELETIONS (NOT DUPLICATIONS)
%==========================================================================
clc; clearvars -except P PGS SVS


SVDATA = SVS.SVDATA;


% Group           Count 
% ----------      ---------------
% <DEL>            110826 (51.0%)   deletion
% <INS:ME:ALU>      44137 (20.3%)   insertion mobile element Alu
% <DUP>             31029 (14.3%)   duplication
% <INS>             15161 (7.0%)    insertion novel sequence
% <INS:ME:LINE1>     7167 (3.3%)    insertion mobile element LINE1
% <INS:ME:SVA>       4154 (1.9%)    insertion mobile element SINE or VNTR or Alu
% <CPX>              3807 (1.8%)    complex structural variant
% <INV>               573 (0.3%)    inversion
% <INS:ME>            447 (0.2%)    insertion mobile element unspecified
% <CTX>                 7 (0.0%)    translocation



REG = 'PROTEIN_CODING__LOF';
TXT = regexp(SVDATA.INFO,REG,'tokens');
CEL = ~cellfun(@isempty, TXT);
ROW = find(CEL);
SVLOF = SVDATA(ROW,:);

% REMOVE SV THAT HAVE LESS THAN 5 INSTANCES
SVLOF(SVLOF.AC < 5 ,:) = [];


% REMOVE SV THAT HAVE OVER 10% HOMALT
SVLOF(SVLOF.FREQ_HOMALT > .10 ,:) = [];




gcounts( SVLOF.ALT )

SVS.SVLOF = SVLOF;



%==========================================================================
%% PULL NUMERICAL DATA OUT OF INFO FIELD
%==========================================================================
clc; clearvars -except P PGS SVS


SVLOF  = SVS.SVLOF;
INFO   = SVLOF.INFO;








NTBL = table();
NTBL.END                 =   1234567890;
NTBL.AFR_AN              =   1000;
NTBL.AFR_AC              =   1000;
NTBL.AFR_N_BI_GENOS      =   1000;
NTBL.AFR_N_HOMREF        =   1000;
NTBL.AFR_N_HET           =   1000;
NTBL.AFR_N_HOMALT        =   1000;
NTBL.AMR_AN              =   1000;
NTBL.AMR_AC              =   1000;
NTBL.AMR_N_BI_GENOS      =   1000;
NTBL.AMR_N_HOMREF        =   1000;
NTBL.AMR_N_HET           =   1000;
NTBL.AMR_N_HOMALT        =   1000;
NTBL.EAS_AN              =   1000;
NTBL.EAS_AC              =   1000;
NTBL.EAS_N_BI_GENOS      =   1000;
NTBL.EAS_N_HOMREF        =   1000;
NTBL.EAS_N_HET           =   1000;
NTBL.EAS_N_HOMALT        =   1000;
NTBL.EUR_AN              =   1000;
NTBL.EUR_AC              =   1000;
NTBL.EUR_N_BI_GENOS      =   1000;
NTBL.EUR_N_HOMREF        =   1000;
NTBL.EUR_N_HET           =   1000;
NTBL.EUR_N_HOMALT        =   1000;
NTBL.OTH_AN              =   1000;
NTBL.OTH_AC              =   1000;
NTBL.OTH_N_BI_GENOS      =   1000;
NTBL.OTH_N_HOMREF        =   1000;
NTBL.OTH_N_HET           =   1000;
NTBL.OTH_N_HOMALT        =   1000;




DTBL = table();
DTBL.AFR_AF              =   "0_01";
DTBL.AFR_FREQ_HOMREF     =   "0_01";
DTBL.AFR_FREQ_HET        =   "0_01";
DTBL.AFR_FREQ_HOMALT     =   "0_01";
DTBL.AMR_AF              =   "0_01";
DTBL.AMR_FREQ_HOMREF     =   "0_01";
DTBL.AMR_FREQ_HET        =   "0_01";
DTBL.AMR_FREQ_HOMALT     =   "0_01";
DTBL.EAS_AF              =   "0_01";
DTBL.EAS_FREQ_HOMREF     =   "0_01";
DTBL.EAS_FREQ_HET        =   "0_01";
DTBL.EAS_FREQ_HOMALT     =   "0_01";
DTBL.EUR_AF              =   "0_01";
DTBL.EUR_FREQ_HOMREF     =   "0_01";
DTBL.EUR_FREQ_HET        =   "0_01";
DTBL.EUR_FREQ_HOMALT     =   "0_01";
DTBL.OTH_AF              =   "0_01";
DTBL.OTH_FREQ_HOMREF     =   "0_01";
DTBL.OTH_FREQ_HET        =   "0_01";
DTBL.OTH_FREQ_HOMALT     =   "0_01";




%%
% THE FOLLOWING REGEXP CAN BE USED TO CAPTURE INTEGERS AND 
% REAL NUMBERS (VALUES WITH DECIMALS), IF ALL DECIMALS ARE FIRST
% REPLACED WITH UNDERSCORES.
% 
% >  SOME_FIELD=(\w+)
% 
% AFTER REPLACING DECIMALS WITH UNDERSCORES, THEY CAN BE CONVERTED
% BACK TO DECIMALS AND THEN CONVERTED FROM CLASS STRING TO NUMERIC
%%
clc; clearvars -except P PGS SVS SVLOF INFO NTBL DTBL

N = numel(INFO);

NTAB = NTBL; NTAB.END(N)    = 10;
DTAB = DTBL; DTAB.AFR_AF(N) = "0_01";




NFIELDS = string(NTAB.Properties.VariableNames');
DFIELDS = string(DTAB.Properties.VariableNames');



for i = 1:N
disp(i)

    ROW = string(INFO(i));


    for j = 1:numel(NFIELDS)

        F   = NFIELDS(j);
        REG = [char(NFIELDS(j)) '=(\w+)'];
        TXT = regexp(ROW,REG,'tokens');

        NTAB.(F)(i) = double(string(TXT));

    end

    WOR = regexprep(ROW,"\.","_");
    
    for j = 1:numel(DFIELDS)

        F   = DFIELDS(j);
        REG = [char(DFIELDS(j)) '=(\w+)'];
        TXT = regexp(WOR,REG,'tokens');

        DTAB.(F)(i) = string(TXT);

    end



end


for j = 1:numel(DFIELDS)

    TXT = regexprep(DTAB.(DFIELDS(j)),"_","\.");
    DTAB.(DFIELDS(j)) = double(TXT);

end




SVDAT = [SVLOF NTAB DTAB];





%==========================================================================
%% PULL OUT GENE NAMES FROM INFO FIELD
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


INFO = SVDAT.INFO;

N = numel(INFO);

GTAB = table();
GTAB.GENE = repmat("NA",N,1);
GTAB.GENZ = repmat("NA",N,1);


for i = 1:N
disp(i)

    ROW = string(INFO(i));
	REG = ['PROTEIN_CODING__LOF=(\w+)'];
	TXT = regexp(ROW,REG,'tokens');
    GTAB.GENE(i) = string(TXT);


    ROW = regexprep(ROW,",","_");
	REG = ['PROTEIN_CODING__LOF=(\w+)'];
	TXT = regexp(ROW,REG,'tokens');
    GTAB.GENZ(i) = string(TXT);

end



SVDAT = [SVDAT GTAB];
%==========================================================================
%% REORGANIZE TABLE VARIABLES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SV = SVDAT;


SV.Properties.VariableNames{'END'} = 'POS_END';
SV.Properties.VariableNames{'ID'} = 'VTYPE';


SV = movevars(SV,{'POS_END'},'After','POS');
SV = movevars(SV,{'VTYPE'},'After','ALT');
SV = movevars(SV,{'GENE'},'After','ALT');
SV = movevars(SV,{'QUAL','FILTER','INFO'},'After','GENZ');


SVDAT = SV;







%==========================================================================
%% COMPUTE ALLELE FREQUENCIES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT
% 
% At single locus with two alleles R and A with frequencies fR and fA,
% respectively, the expected genotype frequencies under random mating are:
% 
% fRR = fR^2
% fAA = fA^2
% fRA = 2*fR*fA
%
% N.B. HARDY-WEINBERG DOES NOT WORK FOR SEX CHROMOSOMES


SV = SVDAT;


% REF COUNT
SV.RC = SV.AN - SV.AC;



% GIVEN...
N   = SV.AN;
A   = SV.AC;
R   = SV.RC;
AA  = SV.N_HOMALT;



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


SV.N   = N;
SV.P   = PPL;
SV.R   = R;
SV.A   = A;
SV.AA  = AA;
SV.RR  = RR;
SV.RA  = RA;
SV.fR  = fR;
SV.fA  = fA;
SV.fRR = fRR;
SV.fAA = fAA;
SV.fRA = fRA;
SV.HW_fRR = HW_fRR;
SV.HW_fAA = HW_fAA;
SV.HW_fRA = HW_fRA;



SVDAT = SV;



%==========================================================================
%% REORGANIZE TABLE VARIABLES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SV = SVDAT;

SV = GENEijk(SV);

SV = movevars(SV,{'RC'},'After','AN');


SV.AN_raw       = SV.AN;
SV.AN_popmax    = SV.AN;
SV.AN_female    = round(SV.AN ./ 2);
SV.AN_male      = round(SV.AN ./ 2);

SV.Properties.VariableNames{'AFR_AN'} = 'AN_afr';
SV.Properties.VariableNames{'AMR_AN'} = 'AN_amr';
SV.Properties.VariableNames{'EAS_AN'} = 'AN_eas';
SV.Properties.VariableNames{'EUR_AN'} = 'AN_nfe';
SV.Properties.VariableNames{'OTH_AN'} = 'AN_oth';

SV.AN_eas_jpn  = SV.AN_eas;
SV.AN_eas_kor  = SV.AN_eas;
SV.AN_eas_oea  = SV.AN_eas;
SV.AN_sas      = SV.AN_eas;
SV.AN_nfe_bgr  = SV.AN_nfe;
SV.AN_nfe_est  = SV.AN_nfe;
SV.AN_nfe_nwe  = SV.AN_nfe;
SV.AN_nfe_onf  = SV.AN_nfe;
SV.AN_nfe_seu  = SV.AN_nfe;
SV.AN_nfe_swe  = SV.AN_nfe;
SV.AN_asj      = SV.AN_nfe;
SV.AN_fin      = SV.AN_nfe;




SV.AC_raw       = SV.AC;
SV.AC_popmax    = SV.AC;
SV.AC_female    = round(SV.AN ./ 2);
SV.AC_male      = round(SV.AN ./ 2);

SV.Properties.VariableNames{'AFR_AC'} = 'AC_afr';
SV.Properties.VariableNames{'AMR_AC'} = 'AC_amr';
SV.Properties.VariableNames{'EAS_AC'} = 'AC_eas';
SV.Properties.VariableNames{'EUR_AC'} = 'AC_nfe';
SV.Properties.VariableNames{'OTH_AC'} = 'AC_oth';


SV.AC_eas_jpn  = SV.AC_eas;
SV.AC_eas_kor  = SV.AC_eas;
SV.AC_eas_oea  = SV.AC_eas;
SV.AC_sas      = SV.AC_eas;
SV.AC_nfe_bgr  = SV.AC_nfe;
SV.AC_nfe_est  = SV.AC_nfe;
SV.AC_nfe_nwe  = SV.AC_nfe;
SV.AC_nfe_onf  = SV.AC_nfe;
SV.AC_nfe_seu  = SV.AC_nfe;
SV.AC_nfe_swe  = SV.AC_nfe;
SV.AC_asj      = SV.AC_nfe;
SV.AC_fin      = SV.AC_nfe;





SV = movevars(SV,{'AN_raw','AN_popmax','AN_female','AN_male'},'After','AF');

SV = movevars(SV,  {'AN_afr','AN_amr','AN_asj','AN_eas','AN_eas_jpn',...
                    'AN_eas_kor','AN_eas_oea','AN_fin','AN_nfe','AN_nfe_bgr',...
                    'AN_nfe_est','AN_nfe_nwe','AN_nfe_onf','AN_nfe_seu',...
                    'AN_nfe_swe','AN_oth','AN_sas'},'After','AN_male');

SV = movevars(SV,{'AC_raw','AC_popmax','AC_female','AC_male'},'After','AN_sas');

SV = movevars(SV,  {'AC_afr','AC_amr','AC_asj','AC_eas','AC_eas_jpn',...
                    'AC_eas_kor','AC_eas_oea','AC_fin','AC_nfe','AC_nfe_bgr',...
                    'AC_nfe_est','AC_nfe_nwe','AC_nfe_onf','AC_nfe_seu',...
                    'AC_nfe_swe','AC_oth','AC_sas'},'After','AC_male');




clc; disp(char(string(SV.Properties.VariableNames')));




SVDAT = SV;



%==========================================================================
%% REORGANIZE TABLE VARIABLES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SV = SVDAT;



SV.AF_raw       = SV.AF;
SV.AF_popmax    = SV.AF;
SV.AF_female    = SV.AF;
SV.AF_male      = SV.AF;

SV.Properties.VariableNames{'N_HOMALT'} = 'nhomalt';

SV.nhomalt_raw      = SV.nhomalt;
SV.nhomalt_popmax   = SV.nhomalt;
SV.nhomalt_female   = SV.nhomalt;
SV.nhomalt_male     = SV.nhomalt;



SV.Properties.VariableNames{'AFR_AF'} = 'AF_afr';
SV.Properties.VariableNames{'AMR_AF'} = 'AF_amr';
SV.Properties.VariableNames{'EAS_AF'} = 'AF_eas';
SV.Properties.VariableNames{'EUR_AF'} = 'AF_nfe';
SV.Properties.VariableNames{'OTH_AF'} = 'AF_oth';

SV.AF_eas_jpn  = SV.AF_eas;
SV.AF_eas_kor  = SV.AF_eas;
SV.AF_eas_oea  = SV.AF_eas;
SV.AF_sas      = SV.AF_eas;
SV.AF_nfe_bgr  = SV.AF_nfe;
SV.AF_nfe_est  = SV.AF_nfe;
SV.AF_nfe_nwe  = SV.AF_nfe;
SV.AF_nfe_onf  = SV.AF_nfe;
SV.AF_nfe_seu  = SV.AF_nfe;
SV.AF_nfe_swe  = SV.AF_nfe;
SV.AF_asj      = SV.AF_nfe;
SV.AF_fin      = SV.AF_nfe;



SV.faf95          = SV.AF;
SV.faf95_afr      = SV.AF_afr;
SV.faf95_amr      = SV.AF_amr;
SV.faf95_eas      = SV.AF_eas;
SV.faf95_nfe      = SV.AF_nfe;
SV.faf95_sas      = SV.AF_sas;
SV.faf99          = SV.AF;
SV.faf99_afr      = SV.AF_afr;
SV.faf99_amr      = SV.AF_amr;
SV.faf99_eas      = SV.AF_eas;
SV.faf99_nfe      = SV.AF_nfe;
SV.faf99_sas      = SV.AF_sas;
SV.n_alt_alleles  = ones(height(SV),1);





SV.Properties.VariableNames{'AFR_N_HOMALT'} = 'nhomalt_afr';
SV.Properties.VariableNames{'AMR_N_HOMALT'} = 'nhomalt_amr';
SV.Properties.VariableNames{'EAS_N_HOMALT'} = 'nhomalt_eas';
SV.Properties.VariableNames{'EUR_N_HOMALT'} = 'nhomalt_nfe';
SV.Properties.VariableNames{'OTH_N_HOMALT'} = 'nhomalt_oth';


SV.nhomalt_eas_jpn  = SV.nhomalt_eas;
SV.nhomalt_eas_kor  = SV.nhomalt_eas;
SV.nhomalt_eas_oea  = SV.nhomalt_eas;
SV.nhomalt_sas      = SV.nhomalt_eas;
SV.nhomalt_nfe_bgr  = SV.nhomalt_nfe;
SV.nhomalt_nfe_est  = SV.nhomalt_nfe;
SV.nhomalt_nfe_nwe  = SV.nhomalt_nfe;
SV.nhomalt_nfe_onf  = SV.nhomalt_nfe;
SV.nhomalt_nfe_seu  = SV.nhomalt_nfe;
SV.nhomalt_nfe_swe  = SV.nhomalt_nfe;
SV.nhomalt_asj      = SV.nhomalt_nfe;
SV.nhomalt_fin      = SV.nhomalt_nfe;



SV.AA_MAF   = SV.AF;
SV.AFR_MAF  = SV.AF_afr;
SV.AMR_MAF  = SV.AF_amr;
SV.EA_MAF   = SV.AF_eas;
SV.EAS_MAF  = SV.AF_eas;
SV.EUR_MAF  = SV.AF_nfe;
SV.SAS_MAF  = SV.AF_eas;
SV.GMAF     = SV.AF;





SV = movevars(SV,{'AF_raw','AF_popmax','AF_female','AF_male'},'After','AC_sas');
SV = movevars(SV,{'nhomalt_raw','nhomalt_popmax','nhomalt_female','nhomalt_male'},'After','AF_male');

SV = movevars(SV,  {'AF_afr','AF_amr','AF_asj','AF_eas','AF_eas_jpn',...
                    'AF_eas_kor','AF_eas_oea','AF_fin','AF_nfe','AF_nfe_bgr',...
                    'AF_nfe_est','AF_nfe_nwe','AF_nfe_onf','AF_nfe_seu',...
                    'AF_nfe_swe','AF_oth','AF_sas'},'After','nhomalt_male');

SV = movevars(SV,  {'faf95','faf95_afr','faf95_amr','faf95_eas','faf95_nfe','faf95_sas',...
                    'faf99','faf99_afr','faf99_amr','faf99_eas','faf99_nfe','faf99_sas',...
                    'n_alt_alleles'},'After','AF_sas');


SV = movevars(SV,  {'nhomalt','nhomalt_afr','nhomalt_amr','nhomalt_asj','nhomalt_eas','nhomalt_eas_jpn',...
                    'nhomalt_eas_kor','nhomalt_eas_oea','nhomalt_fin','nhomalt_nfe','nhomalt_nfe_bgr',...
                    'nhomalt_nfe_est','nhomalt_nfe_nwe','nhomalt_nfe_onf','nhomalt_nfe_seu',...
                    'nhomalt_nfe_swe','nhomalt_oth','nhomalt_sas'},'After','n_alt_alleles');


SV = movevars(SV,  {'AA_MAF','AFR_MAF','AMR_MAF','EA_MAF',...
                    'EAS_MAF','EUR_MAF','SAS_MAF','GMAF'},'After','nhomalt_sas');


clc; disp(char(string(SV.Properties.VariableNames')));




SVDAT = SV;





%==========================================================================
%% REORGANIZE TABLE VARIABLES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SV = SVDAT;


N = height(SV);
ONES = ones(N,1);
NANS = nan(N,1);
NaNs = repmat("NaN",N,1);
NA   = repmat("NA",N,1);
MISS = repmat(PGS.GNOMAD.EXON(1),N,1);


SV.pab_max                  = ONES;
SV.popmax                   = NA;
SV.ALLELE_NUM               = ONES;
SV.APPRIS                   = NaNs;
SV.BIOTYPE                  = repmat("protein_coding",N,1);
SV.CANONICAL                = ONES;
SV.CCDS                     = NaNs;
SV.CLIN_SIG                 = NaNs;
SV.Consequence              = SV.INFO;
SV.DISTANCE                 = NaNs;
SV.ENSP                     = NaNs;
SV.Existing_variation       = NaNs;
SV.EXON                     = MISS;
SV.Feature                  = NaNs;
SV.Feature_type             = repmat("Transcript",N,1);
SV.FLAGS                    = NaNs;
SV.GENE_ID                  = SV.GENZ;
SV.GENE_PHENO               = ONES-1;
SV.HGNC_ID                  = NaNs;
SV.HIGH_INF_POS             = NaNs;
SV.IMPACT                   = repmat("High",N,1);
SV.INTRON                   = NaNs;
SV.MINIMISED                = ONES;
SV.MOTIF_NAME               = NaNs;
SV.MOTIF_POS                = NaNs;
SV.MOTIF_SCORE_CHANGE       = NaNs;
SV.PHENO                    = NaNs;
SV.PolyPhen                 = NaNs;
SV.Protein_position         = NaNs;
SV.PUBMED                   = NaNs;
SV.SIFT                     = NaNs;
SV.SOMATIC                  = NaNs;
SV.STRAND                   = ONES;
SV.SWISSPROT                = SV.GENZ;
SV.SYMBOL_SOURCE            = repmat("HGNC",N,1);
SV.TREMBL                   = NaNs;
SV.TSL                      = NaNs;
SV.UNIPARC                  = NaNs;
SV.VARIANT_CLASS            = repmat("SV",N,1);
SV.BaseQRankSum             = NANS;
SV.ClippingRankSum          = NANS;
SV.DP                       = NANS;
SV.filter                   = SV.FILTER;
SV.FS                       = NANS;
SV.has_star                 = NANS;
SV.InbreedingCoeff          = NANS;
SV.MQ                       = ONES+59;
SV.MQRankSum                = ONES-.5;
SV.names                    = NA;
SV.nonpar                   = NANS;
SV.QD                       = NANS;
SV.quality                  = SV.QUAL;
SV.ReadPosRankSum           = NANS;
SV.rf_label                 = NaNs;
SV.rf_negative_label        = NANS;
SV.rf_positive_label        = NANS;
SV.rf_tp_probability        = NANS;
SV.rf_train                 = NANS;
SV.SOR                      = NANS;
SV.transmitted_singleton    = NaNs;
SV.VQSLOD                   = NaNs;
SV.VQSR_culprit             = NA;
SV.VQSR_NEGATIVE_TRAIN_SITE = NANS;
SV.VQSR_POSITIVE_TRAIN_SITE = MISS;
SV.was_mixed                = NANS;
SV.allele_type              = repmat("SV",N,1);
SV.LoF                      = repmat("HC",N,1);
SV.LoF_filter               = NaNs;
SV.LoF_flags                = NaNs;
SV.segdup                   = NANS;
SV.lcr                      = NANS;
SV.decoy                    = NANS;


SV.FILTER = [];
SV.INFO   = [];
SV.QUAL   = [];
SV.GENZ   = [];


SV = movevars(SV,  {'N','P','R','A','AA','RR',...
                    'RA','fR','fA','fRR','fAA','fRA',...
                    'HW_fRR','HW_fAA','HW_fRA'},'After','decoy');


clc; disp(char(string(SV.Properties.VariableNames')));


SVDAT = SV;



%==========================================================================
%% REORGANIZE TABLE VARIABLES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SV = SVDAT;



N = height(SV);
NANS = nan(N,1);



SV.VCR              = nan(N,1);
SV.VCR_raw          = nan(N,1);
SV.VCR_popmax       = nan(N,1);
SV.VCR_afr          = nan(N,1);
SV.VCR_amr          = nan(N,1);
SV.VCR_asj          = nan(N,1);
SV.VCR_eas          = nan(N,1);
SV.VCR_eas_jpn      = nan(N,1);
SV.VCR_eas_kor      = nan(N,1);
SV.VCR_eas_oea      = nan(N,1);
SV.VCR_fin          = nan(N,1);
SV.VCR_nfe          = nan(N,1);
SV.VCR_nfe_bgr      = nan(N,1);
SV.VCR_nfe_est      = nan(N,1);
SV.VCR_nfe_nwe      = nan(N,1);
SV.VCR_nfe_onf      = nan(N,1);
SV.VCR_nfe_seu      = nan(N,1);
SV.VCR_nfe_swe      = nan(N,1);
SV.VCR_oth          = nan(N,1);
SV.VCR_sas          = nan(N,1);
SV.GCR              = nan(N,1);
SV.GCR_raw          = nan(N,1);
SV.GCR_popmax       = nan(N,1);
SV.GCR_afr          = nan(N,1);
SV.GCR_amr          = nan(N,1);
SV.GCR_asj          = nan(N,1);
SV.GCR_eas          = nan(N,1);
SV.GCR_eas_jpn      = nan(N,1);
SV.GCR_eas_kor      = nan(N,1);
SV.GCR_eas_oea      = nan(N,1);
SV.GCR_fin          = nan(N,1);
SV.GCR_nfe          = nan(N,1);
SV.GCR_nfe_bgr      = nan(N,1);
SV.GCR_nfe_est      = nan(N,1);
SV.GCR_nfe_nwe      = nan(N,1);
SV.GCR_nfe_onf      = nan(N,1);
SV.GCR_nfe_seu      = nan(N,1);
SV.GCR_nfe_swe      = nan(N,1);
SV.GCR_oth          = nan(N,1);
SV.GCR_sas          = nan(N,1);
SV.isLOF            = ones(N,1);





SV = movevars(SV,  {'N','P','R','A','AA','RR',...
                    'RA','fR','fA','fRR','fAA','fRA',...
                    'HW_fRR','HW_fAA','HW_fRA'},'After','decoy');


clc; disp(char(string(SV.Properties.VariableNames')));


SVDAT = SV;







%==========================================================================
%% SAVE A FULL VERSION OF THE SVDAT
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SVS.SVFULL = SVDAT;






%==========================================================================
%% REORGANIZE TABLE VARIABLES
%==========================================================================
clc; clearvars -except P PGS SVS SVDAT


SV = SVDAT;


SV.N_BI_GENOS 		= [];
SV.N_HOMREF 		= [];
SV.N_HET 			= [];
SV.FREQ_HOMREF 		= [];
SV.FREQ_HET 		= [];
SV.FREQ_HOMALT 		= [];
SV.AFR_N_BI_GENOS 	= [];
SV.AFR_N_HOMREF 	= [];
SV.AFR_N_HET 		= [];
SV.AMR_N_BI_GENOS 	= [];
SV.AMR_N_HOMREF 	= [];
SV.AMR_N_HET 		= [];
SV.EAS_N_BI_GENOS 	= [];
SV.EAS_N_HOMREF 	= [];
SV.EAS_N_HET 		= [];
SV.EUR_N_BI_GENOS 	= [];
SV.EUR_N_HOMREF 	= [];
SV.EUR_N_HET 		= [];
SV.OTH_N_BI_GENOS 	= [];
SV.OTH_N_HOMREF 	= [];
SV.OTH_N_HET 		= [];
SV.AFR_FREQ_HOMREF 	= [];
SV.AFR_FREQ_HET 	= [];
SV.AFR_FREQ_HOMALT 	= [];
SV.AMR_FREQ_HOMREF 	= [];
SV.AMR_FREQ_HET 	= [];
SV.AMR_FREQ_HOMALT 	= [];
SV.EAS_FREQ_HOMREF 	= [];
SV.EAS_FREQ_HET 	= [];
SV.EAS_FREQ_HOMALT 	= [];
SV.EUR_FREQ_HOMREF 	= [];
SV.EUR_FREQ_HET 	= [];
SV.EUR_FREQ_HOMALT 	= [];
SV.OTH_FREQ_HOMREF 	= [];
SV.OTH_FREQ_HET 	= [];
SV.OTH_FREQ_HOMALT	= [];


clc; disp(char(string(SV.Properties.VariableNames')));




SVLOF = SV;

SVS.SVLOF = SVLOF;

%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS SVS


DATAS  = PGS.DATAS;
DAWES  = PGS.DAWES;
GNOMAD = PGS.GNOMAD;
SVFULL = SVS.SVFULL;
SVLOF  = SVS.SVLOF;

save([P.mat P.f 'PGS_STEP052_OUTPUT.mat'],'GNOMAD','DATAS','DAWES','SVFULL','SVLOF');

clc; clearvars -except P PGS SVS

