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



PGS = load([P.mat P.f 'PGS_STEP046_OUTPUT.mat']);


%==========================================================================
%% IMPORT DAWES DATA
%==========================================================================
clc; clearvars -except P PGS



% HGNC = readtable([P.csv P.f 'HGNC_RESULTS_FOR_DAWES_GENES.csv']);
% 
%--------------------------------------------------------
% DAWES_S2.xlsx
%   S0_Legend
%   S1_GD_informatics_toolkit
%   S2_cell_ess_mouse_nonlethal
%   S3_candidate_lethal_genes
%   S4_extra_mouselethalgenes




DAWES2_S0 = readtable([P.csv '/dawes/DAWES_DATA.xlsx'], 'Sheet','LEGEND');
DAWES2_S1 = readtable([P.csv '/dawes/DAWES_DATA.xlsx'], 'Sheet','S1_GD_informatics_toolkit');
DAWES2_S2 = readtable([P.csv '/dawes/DAWES_DATA.xlsx'], 'Sheet','S2_cell_ess_mouse_nonlethal');
DAWES2_S3 = readtable([P.csv '/dawes/DAWES_DATA.xlsx'], 'Sheet','S3_candidate_lethal_genes');
DAWES2_S4 = readtable([P.csv '/dawes/DAWES_DATA.xlsx'], 'Sheet','S4_extra_mouselethalgenes');

DAWES2_S0 = convertvars(DAWES2_S0,@iscell,'string');
DAWES2_S1 = convertvars(DAWES2_S1,@iscell,'string');
DAWES2_S2 = convertvars(DAWES2_S2,@iscell,'string');
DAWES2_S3 = convertvars(DAWES2_S3,@iscell,'string');
DAWES2_S4 = convertvars(DAWES2_S4,@iscell,'string');




%--------------------------------------------------------
% DAWES_S3.xlsx
%   S0_Legend
%   S1_lethal_search_terms
%   S2_lethal_search_hit_details
%   S3_cell_KO_tallies
%   S4_lethal_MP_terms

DAWES3_S0 = readtable([P.csv '/dawes/DAWES_SUPP.xlsx'], 'Sheet','S0_Legend');
DAWES3_S1 = readtable([P.csv '/dawes/DAWES_SUPP.xlsx'], 'Sheet','S1_lethal_search_terms');
DAWES3_S2 = readtable([P.csv '/dawes/DAWES_SUPP.xlsx'], 'Sheet','S2_lethal_search_hit_details');
DAWES3_S3 = readtable([P.csv '/dawes/DAWES_SUPP.xlsx'], 'Sheet','S3_cell_KO_tallies');
DAWES3_S4 = readtable([P.csv '/dawes/DAWES_SUPP.xlsx'], 'Sheet','S4_lethal_MP_terms');

DAWES3_S0 = convertvars(DAWES3_S0,@iscell,'string');
DAWES3_S1 = convertvars(DAWES3_S1,@iscell,'string');
DAWES3_S2 = convertvars(DAWES3_S2,@iscell,'string');
DAWES3_S3 = convertvars(DAWES3_S3,@iscell,'string');
DAWES3_S4 = convertvars(DAWES3_S4,@iscell,'string');



%==========================================================================
%% SAVE DAWES DATA
%==========================================================================
clc; clearvars -except P PGS...
DAWES2_S0 DAWES2_S1 DAWES2_S2 DAWES2_S3 DAWES2_S4 ...
DAWES3_S0 DAWES3_S1 DAWES3_S2 DAWES3_S3 DAWES3_S4




save([P.mat P.f 'DAWES_S2_S3.mat'],...
    'DAWES2_S0','DAWES2_S1','DAWES2_S2','DAWES2_S3','DAWES2_S4',...
    'DAWES3_S0','DAWES3_S1','DAWES3_S2','DAWES3_S3','DAWES3_S4');



DAWES.DAWES2_S0 = DAWES2_S0;
DAWES.DAWES2_S1 = DAWES2_S1;
DAWES.DAWES2_S2 = DAWES2_S2;
DAWES.DAWES2_S3 = DAWES2_S3;
DAWES.DAWES2_S4 = DAWES2_S4;

DAWES.DAWES3_S0 = DAWES3_S0;
DAWES.DAWES3_S1 = DAWES3_S1;
DAWES.DAWES3_S2 = DAWES3_S2;
DAWES.DAWES3_S3 = DAWES3_S3;
DAWES.DAWES3_S4 = DAWES3_S4;


DATAS = PGS.DATAS;
GNOMAD = PGS.GNOMAD;

save([P.mat P.f 'PGS_STEP048_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');










