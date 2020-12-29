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
%------------------------------------------------










%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP055_B_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;   % GNOMAD DATASET




% GNOMAD COLUMNS   1-111
% PLI    COLUMNS 112-121
% DAWES  COLUMNS 122-183
% MK1    COLUMNS 184-216
% MK2    COLUMNS 217-276
% SCR    COLUMNS 277-279
% TARG   COLUMNS 280-286
%==========================================================================
%% REMOVE SITES WITH FEW SEQUENCED ALLELES (AN)
%==========================================================================
clc; clearvars -except P PGS GNOMAD



N = round(quantile(GNOMAD.AN,.01));

GNOMAD( GNOMAD.AN < N , :) = [];





Nreport(GNOMAD)
%==========================================================================
%% REMOVE SEX CHROMOSOMES
%==========================================================================
clc; clearvars -except P PGS GNOMAD



GNOMAD(GNOMAD.CHR == 23 , :) = [];









Nreport(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS THAT HAVE OVER 50% MAF
%==========================================================================
clc; clearvars -except P PGS GNOMAD





GNOMAD(GNOMAD.AF > .50 , :) = [];






Nreport(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS THAT HAVE MORE THAN N% HOMOZYGOUS MINOR ALLELES
%==========================================================================
clc; clearvars -except P PGS GNOMAD




N = .05;

i = (GNOMAD.NHOMALT ./ GNOMAD.AC) > N;

GNOMAD(i , :) = [];







Nreport(GNOMAD)
%==========================================================================
%% REMOVE ROWS FLAGGED BY GNOMAD PIPELINE
%==========================================================================
clc; clearvars -except P PGS GNOMAD





gcounts( GNOMAD.inSegDup )
gcounts( GNOMAD.inLowComplex )
gcounts( GNOMAD.inDecoy )
gcounts( GNOMAD.lof_flags ~= "NaN" )




GNOMAD(GNOMAD.inSegDup == 1 , :)        = [];
GNOMAD(GNOMAD.inLowComplex == 1 , :)    = [];
GNOMAD(GNOMAD.inDecoy == 1 , :)         = [];
GNOMAD(GNOMAD.lof_flags ~= "NaN" , :)   = [];





Nreport(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS WITH MAF > N
%==========================================================================
clc; clearvars -except P PGS GNOMAD




%N = .20;

N = .05;

GNOMAD(GNOMAD.AF > N ,:) = [];








Nreport(GNOMAD)
%==========================================================================
%% ONLY USE CANONICAL TRANSCRIPTS
%==========================================================================
% clc; clearvars -except P PGS GNOMAD
% 
% 
% GNOMAD( GNOMAD.canonical ~= 1, : ) = [];
% 
% DMAD( DMAD.canonical ~= 1, : ) = [];
% 
% 
% height(GNOMAD)
% height(DMAD)
% 
% 
%==========================================================================
%% REMOVE SITES THAT APPEAR LIKE DUPLICATE ALLELES (RUN AS LAST EXCLUSION)
%==========================================================================
clc; clearvars -except P PGS GNOMAD




GNOM = GNOMAD;
[C,ia,ic] = unique([GNOM.CHR GNOM.POS GNOM.AC],'rows','stable');
GNOM = GNOM(ia,:);
GNOMAD = GNOM;


GNOMAD = GENEijk(GNOMAD);






Nreport(GNOMAD)
%==========================================================================
%% IMPORT OMIM DATA
%==========================================================================
clc; clearvars -except P PGS GNOMAD



% ABOUT OMIM
%---------------------------------------------
%{
Welcome to OMIM®, Online Mendelian Inheritance in Man®. OMIM is a
comprehensive, authoritative compendium of human genes and genetic
phenotypes that is freely available and updated daily. The full-text,
referenced overviews in OMIM contain information on all known mendelian
disorders and over 15,000 genes. OMIM focuses on the relationship between
phenotype and genotype. It is updated daily, and the entries contain
copious links to other genetics resources.

This database was initiated in the early 1960s by Dr. Victor A. McKusick as
a catalog of mendelian traits and disorders, entitled Mendelian Inheritance
in Man (MIM). Twelve book editions of MIM were published between 1966 and
1998. The online version, OMIM, was created in 1985 by a collaboration
between the National Library of Medicine and the William H. Welch Medical
Library at Johns Hopkins. It was made generally available on the internet
starting in 1987. In 1995, OMIM was developed for the World Wide Web by
NCBI, the National Center for Biotechnology Information.

OMIM is authored and edited at the McKusick-Nathans Institute of Genetic
Medicine, Johns Hopkins University School of Medicine, under the direction
of Dr. Ada Hamosh.


https://www.omim.org/about

https://www.omim.org/static/omim/data/mim2gene.txt

%}
%---------------------------------------------



P.mim = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2/PGS_CSV/omim/mim2gene.txt';

OPT = detectImportOptions(P.mim);
TBL = readtable(P.mim,OPT);


[ai,bi] = ismember(GNOMAD.GENE,TBL.GeneSym);


GNOMAD.MIMNumber = nan(height(GNOMAD),1);
GNOMAD.MIMNumber(ai) = TBL.MIMNumber(bi(ai),:);


sum(~isnan(GNOMAD.MIMNumber)) / height(GNOMAD)

sum(GNOMAD.HAS_OMIM) / height(GNOMAD)



%return
%==========================================================================
%% IMPORT CLINVAR DATA
%==========================================================================
clc; clearvars -except P PGS GNOMAD




% ABOUT CLINVAR
%---------------------------------------------
%{
ClinVar is a a freely accessible, public archive of reports of the
relationships among human variations and phenotypes hosted by the National
Center for Biotechnology Information (NCBI) and funded by intramural
National Institutes of Health (NIH) funding.


https://www.clinicalgenome.org/data-sharing/clinvar/

https://ftp.ncbi.nlm.nih.gov/pub/clinvar/

%}
%---------------------------------------------


P.clin = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2/PGS_CSV/omim/clinvar.txt';

OPT = detectImportOptions(P.clin);
TBL = readtable(P.clin,OPT);

TBL = convertvars(TBL,@iscell,'string');


GNOMAD.CHRPOS = uint64(GNOMAD.CHR .* 10000000000 + GNOMAD.POS);
TBL.CHRPOS    = uint64(TBL.chrom .* 10000000000 + TBL.pos);


[ai,bi] = ismember(GNOMAD.CHRPOS,TBL.CHRPOS);


GNOMAD.CLINVAR_PATHOGENIC           = nan(height(GNOMAD),1);
GNOMAD.CLINVAR_LIKELY_PATHOGENIC    = nan(height(GNOMAD),1);
GNOMAD.CLINVAR_UNCERTAIN            = nan(height(GNOMAD),1);
GNOMAD.CLINVAR_LIKELY_BENIGN        = nan(height(GNOMAD),1);
GNOMAD.CLINVAR_BENIGN               = nan(height(GNOMAD),1);


GNOMAD.CLINVAR_PATHOGENIC(ai)           = TBL.pathogenic(bi(ai),:);
GNOMAD.CLINVAR_LIKELY_PATHOGENIC(ai)    = TBL.likely_pathogenic(bi(ai),:);
GNOMAD.CLINVAR_UNCERTAIN(ai)            = TBL.uncertain_significance(bi(ai),:);
GNOMAD.CLINVAR_LIKELY_BENIGN(ai)        = TBL.likely_benign(bi(ai),:);
GNOMAD.CLINVAR_BENIGN(ai)               = TBL.benign(bi(ai),:);


sum(~isnan(GNOMAD.CLINVAR_PATHOGENIC)) / height(GNOMAD)



%==========================================================================
%% IMPORT CADD SCORES
%==========================================================================
clc; clearvars -except P PGS GNOMAD



% ABOUT CADD
%---------------------------------------------
%{

Use STEP057_MAKE_VCF.m to create a vcf file, then upload to:

https://cadd.gs.washington.edu/score

Then import the tsv file below.

%}
%---------------------------------------------


P.cadd = '/Volumes/T7/UCSF/PGS/CADD_SCORES.txt';

OPT = detectImportOptions(P.cadd);
CADD = readtable(P.cadd,OPT);


CADD = convertvars(CADD,@iscell,'string');

CADD.CHRPOS = uint64(CADD.x_Chrom .* 10000000000 + CADD.Pos);


[ai,bi] = ismember(GNOMAD.CHRPOS,CADD.CHRPOS);


GNOMAD.CADD_RAW         = nan(height(GNOMAD),1);
GNOMAD.CADD_PHRED       = nan(height(GNOMAD),1);


GNOMAD.CADD_RAW(ai)     = CADD.RawScore(bi(ai),:);
GNOMAD.CADD_PHRED(ai)	= CADD.PHRED(bi(ai),:);


sum(~isnan(GNOMAD.CADD_RAW)) / height(GNOMAD)




%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


GAD = PGS.GNOMAD;

TABLES = PGS.TABLES;

save([P.mat P.f 'PGS_STEP056_OUTPUT.mat'],'GNOMAD','TABLES','GAD');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP056_OUTPUT.csv'])


