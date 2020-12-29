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
%% MAKE VCF FILE FOR CADD SCORING
%==========================================================================
clc; clearvars -except P PGS GNOMAD

% THE VCF FILE SHOULD BE FORMATTED LIKE...
% 
% #CHROM POS ID REF ALT
% 1 1417517 rs772346770 G C
% 1 1417961 rs761985087 C CT
% 1 1418426 rs991754583 G C
% 1 1420548 rs913743808 CCAGGT C

% UPLOAD VCF FILE TO: https://cadd.gs.washington.edu/score



VCF = GNOMAD(:,[4,5,17,18,20]);

VCF.Properties.VariableNames = {'CHROM','POS','ID','REF','ALT'};













