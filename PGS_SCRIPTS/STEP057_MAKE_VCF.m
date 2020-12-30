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
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P

PGS = load([P.mat P.f 'PGS_STEP056_OUTPUT.mat']);



GNOMAD  = PGS.GNOMAD;   % GNOMAD DATASET



peek(GNOMAD)
%==========================================================================
%% PREPARE TABLE USING VCF-STYLE COLUMNS FOR UPLOAD TO CADD ANNOTATION SITE
%==========================================================================
% 
% 
% THE VCF FILE SHOULD BE FORMATTED LIKE...
% 
% #CHROM POS ID REF ALT
% 1 1417517 rs772346770 G C
% 1 1417961 rs761985087 C CT
% 1 1418426 rs991754583 G C
% 1 1420548 rs913743808 CCAGGT C
% 
% 
% UPLOAD VCF FILE TO: 
% 
%       https://cadd.gs.washington.edu/score
% 
% 
%---------------------------------------------
clc; clearvars -except P PGS GNOMAD



VCF = table();

VCF.CHROM = GNOMAD.CHR;
VCF.POS   = GNOMAD.POS;
VCF.ID    = GNOMAD.rs;
VCF.REF   = GNOMAD.refBase;
VCF.ALT   = GNOMAD.altBase;


% clc; numel(unique(GNOMAD.altBase))
% clc; numel(unique(GNOMAD.altAllele))
% 
% VCF = GNOMAD(:,[5,6,18,19,21]);
% VCF.Properties.VariableNames = {'CHROM','POS','ID','REF','ALT'};


VCF.POS(ismissing(VCF.POS)) = ".";
VCF.ID(ismissing(VCF.ID))   = ".";
VCF.REF(ismissing(VCF.REF)) = ".";
VCF.ALT(ismissing(VCF.ALT)) = ".";

VCF.REF(VCF.REF == "0") = ".";
VCF.ALT(VCF.ALT == "0") = ".";
VCF.ID(VCF.ID == "NaN") = ".";


%==========================================================================
%% OUTPUT A SPACE-DELIMITED TEXT FILE
%==========================================================================
% 
% NOTE:
% BEFORE UPLOAD, THIS TEXT FILE MUST BE OPENED, AND A HASHTAG PLACED 
% AT THE START OF THE HEADER LINE, TO READ:
% 
%       #CHROM POS ID REF ALT
% 
% THEN THE FILE EXTENSION SHOULD BE CHANGED FROM .TXT TO .VCF
% AND UPLOADED TO https://cadd.gs.washington.edu/score
% 
%---------------------------------------------
clc; clearvars -except P PGS GNOMAD VCF



writetable(VCF,[P.csv P.f 'cadd/GNOMAD_CADD_VCF_FILE.txt'],'Delimiter',' ')






