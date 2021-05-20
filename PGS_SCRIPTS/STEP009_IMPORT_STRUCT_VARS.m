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

PGS = load([P.mat P.f 'PGS_STEP050_OUTPUT.mat']);



GNOMAD  = PGS.GNOMAD;   % GNOMAD DATASET



peek(GNOMAD)
%==========================================================================
%% IMPORT STRUCTURAL VARIANTS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



% ABOUT SVARS
%---------------------------------------------
%{

Use STEP057_MAKE_VCF.m to create a vcf file, then upload to:

https://cadd.gs.washington.edu/score

Then import the tsv file below.

%}
%---------------------------------------------


P.svars = '/Volumes/T7/UCSF/PGS/STRUCTURAL_VARIANTS.txt';

OPT = detectImportOptions(P.svars);
SVARS = readtable(P.svars,OPT);


SVARS.FILTER = string(SVARS.FILTER);

SVARS(SVARS.FILTER ~= "PASS" ,:) = [];
SVARS(SVARS.QUAL < 500 ,:) = [];

CEL = cellfun(@split, SVARS.INFO,  repmat({[";"]}, size(SVARS.INFO,1) ,1)  , 'UniformOutput',false );



N = size(CEL,1);


for i = 1:N

    CEL{i} = CEL{i}(1:40);

end






%==========================================================================
%% IMPORT STRUCTURAL VARIANTS
%==========================================================================
clc; clearvars -except P PGS GNOMAD SVARS CEL


% 'AN=21694'
% 'AC=2'
% 'AF=9.2e-05'
% 'N_BI_GENOS=10847'
% 'N_HOMREF=10845'
% 'N_HET=2'
% 'N_HOMALT=0'
% 'FREQ_HOMREF=0.999816'
% 'FREQ_HET=0.000184383'
% 'FREQ_HOMALT=0'




N = size(CEL,1);


for i = 1:N

    CEL{i} = string(CEL{i});

end



AN          = nan(N,1);
AC          = nan(N,1);
AF          = nan(N,1);
N_BI_GENOS  = nan(N,1);
N_HOMREF    = nan(N,1);
N_HET       = nan(N,1);
N_HOMALT    = nan(N,1);
FREQ_HOMREF = nan(N,1);
FREQ_HET    = nan(N,1);
FREQ_HOMALT = nan(N,1);


for i = 1:N

    C = CEL{i};

    [~,AN_i]            = regexp(C,'AN=');
    [~,AC_i]            = regexp(C,'AC=');
    [~,AF_i]            = regexp(C,'AF=');
    [~,N_BI_GENOS_i]    = regexp(C,'N_BI_GENOS=');
    [~,N_HOMREF_i]      = regexp(C,'N_HOMREF=');
    [~,N_HET_i]         = regexp(C,'N_HET=');
    [~,N_HOMALT_i]      = regexp(C,'N_HOMALT=');
    [~,FREQ_HOMREF_i]   = regexp(C,'FREQ_HOMREF=');
    [~,FREQ_HET_i]      = regexp(C,'FREQ_HET=');
    [~,FREQ_HOMALT_i]   = regexp(C,'FREQ_HOMALT=');


    for j = 1:40

        if AN_i{j} == 3
            AN(i) = str2double(C{j}(4:end));
        end

        if AC_i{j} == 3
            AC(i) = str2double(C{j}(4:end));
        end

        if AF_i{j} == 3
            AF(i) = str2double(C{j}(4:end));
        end

        if N_BI_GENOS_i{j} == 11
            N_BI_GENOS(i) = str2double(C{j}(12:end));
        end

        if N_HOMREF_i{j} == 9
            N_HOMREF(i) = str2double(C{j}(10:end));
        end

        if N_HET_i{j} == 6
            N_HET(i) = str2double(C{j}(7:end));
        end

        if N_HOMALT_i{j} == 9
            N_HOMALT(i) = str2double(C{j}(10:end));
        end


        if FREQ_HOMREF_i{j} == 12
            FREQ_HOMREF(i) = str2double(C{j}(13:end));
        end


        if FREQ_HET_i{j} == 9
            FREQ_HET(i) = str2double(C{j}(10:end));
        end


        if FREQ_HOMALT_i{j} == 12
            FREQ_HOMALT(i) = str2double(C{j}(13:end));
        end

    end


end



STRUCTVARS = SVARS;

STRUCTVARS.AN = AN;
STRUCTVARS.AC = AC;
STRUCTVARS.AF = AF;
STRUCTVARS.N_BI_GENOS = N_BI_GENOS;
STRUCTVARS.N_HOMREF = N_HOMREF;
STRUCTVARS.N_HET = N_HET;
STRUCTVARS.N_HOMALT = N_HOMALT;
STRUCTVARS.FREQ_HOMREF = FREQ_HOMREF;
STRUCTVARS.FREQ_HET = FREQ_HET;
STRUCTVARS.FREQ_HOMALT = FREQ_HOMALT;


%==========================================================================
%% SAVE DATA
%==========================================================================
clc; clearvars -except P PGS GNOMAD SVARS CEL STRUCTVARS


save([P.mat P.f 'STRUCTVARS.mat'],'STRUCTVARS','SVARS');


