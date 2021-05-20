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


PGS = load([P.mat P.f 'PGS_STEP090_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;






peek(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% IMPORT OMIM DATA
%==========================================================================
% 
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
% 
%---------------------------------------------
clc; clearvars -except P PGS GNOMAD


P.mim = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2/PGS_CSV/omim/mim2gene.txt';

OPT = detectImportOptions(P.mim);
TBL = readtable(P.mim,OPT);


[ai,bi] = ismember(GNOMAD.GENE,TBL.GeneSym);


GNOMAD.OMIM     = nan(height(GNOMAD),1);
GNOMAD.OMIM(ai) = TBL.MIMNumber(bi(ai),:);






%==========================================================================
%% REPORT STATS ON   OMIM - GNOMAD - DAWES   MATCHES
%==========================================================================
clc; clearvars -except P PGS GNOMAD



NUM = sum(~isnan(GNOMAD.OMIM));
PCT = NUM / height(GNOMAD)*100;
fprintf('OMIM ENTRIES MATCHING GNOMAD SNVs:   %7.0f   (%2.1f %%) \n',NUM,PCT)


NUM = sum(~isnan(GNOMAD.DAW_mim));
PCT = NUM / height(GNOMAD)*100;
fprintf('OMIM ENTRIES MATCHING DAWES OMIMS:   %7.0f   (%2.1f %%) \n',NUM,PCT)






%return
%==========================================================================
%% IMPORT CLINVAR DATA
%==========================================================================
% 
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
clc; clearvars -except P PGS GNOMAD




P.clin = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2/PGS_CSV/clinvar/clinvar.txt';

OPT = detectImportOptions(P.clin);
TBL = readtable(P.clin,OPT);

TBL = convertvars(TBL,@iscell,'string');


TBL.CHRPOS = makeCHRPOS(TBL.chrom,TBL.pos);
TBL = movevars(TBL,{'CHRPOS'},'Before','chrom');




% CHECK ALIGNMENT OF GNOMAD.CHRPOS & CLINVAR.CHRPOS
%---------------------------------------------
% LOF = GNOMAD(GNOMAD.isLOF == 1 ,:);
% STR = GNOMAD(GNOMAD.isLOF ~= 1 ,:);
% [Lai,Lbi] = ismember(LOF.CHRPOS+0,TBL.CHRPOS); sum(Lai)
% [Lai,Lbi] = ismember(LOF.CHRPOS+1,TBL.CHRPOS); sum(Lai)
% [Lai,Lbi] = ismember(LOF.CHRPOS-1,TBL.CHRPOS); sum(Lai)
% [Lai,Lbi] = ismember(LOF.CHRPOS+2,TBL.CHRPOS); sum(Lai)
% [Lai,Lbi] = ismember(LOF.CHRPOS-2,TBL.CHRPOS); sum(Lai)
% [Sai,Sbi] = ismember(STR.CHRPOS+0,TBL.CHRPOS); sum(Sai)
% [Sai,Sbi] = ismember(STR.CHRPOS+1,TBL.CHRPOS); sum(Sai)
% [Sai,Sbi] = ismember(STR.CHRPOS-1,TBL.CHRPOS); sum(Sai)
% [Sai,Sbi] = ismember(STR.CHRPOS+2,TBL.CHRPOS); sum(Sai)
% [Sai,Sbi] = ismember(STR.CHRPOS-2,TBL.CHRPOS); sum(Sai)
%---------------------------------------------



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


NUM = sum(~isnan(GNOMAD.CLINVAR_PATHOGENIC));
PCT = NUM / height(GNOMAD)*100;
fprintf('PATHOGENIC CLINVAR SVNs:   %10.0f   (%2.1f %%) \n',NUM,PCT)






% return
%==========================================================================
%% IMPORT CADD SCORES
%==========================================================================
% 
% ABOUT CADD
%---------------------------------------------
%{

Use STEP057_MAKE_VCF.m to create a vcf file, then upload to:

https://cadd.gs.washington.edu/score

Then import the tsv file below.

%}
%---------------------------------------------
clc; clearvars -except P PGS GNOMAD




P.cadd = [P.csv '/cadd/GNOMAD_CADD_OUTPUT_V02.txt'];

OPT = detectImportOptions(P.cadd);
CADD = readtable(P.cadd,OPT);


CADD = convertvars(CADD,@iscell,'string');


CADD.CHRPOS = makeCHRPOS( CADD.Chrom , CADD.Pos );
CADD = movevars(CADD,{'CHRPOS'},'Before','Chrom');




[ai,bi] = ismember(GNOMAD.CHRPOS,CADD.CHRPOS);


GNOMAD.CADD_RAW         = nan(height(GNOMAD),1);
GNOMAD.CADD_PHRED       = nan(height(GNOMAD),1);


GNOMAD.CADD_RAW(ai)     = CADD.RawScore(bi(ai),:);
GNOMAD.CADD_PHRED(ai)	= CADD.PHRED(bi(ai),:);





NUM = sum(~isnan(GNOMAD.CADD_RAW));
PCT = NUM / height(GNOMAD)*100;
fprintf('SVNs WITH CADD ENTRY:   %10.0f   (%2.1f %%) \n',NUM,PCT)






%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DATAS = PGS.DATAS;
DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP100_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');



