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
P.funs   = [P.home filesep 'PGS_FUNS'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------




%==========================================================================
%% GET FULL PATHS TO MAT FILES
%==========================================================================
clearvars -except P


WFILES = dir([P.gnomat P.f '**/*.mat']);

FILES = strcat(string({WFILES.folder})',...
            repmat(string(filesep),numel(WFILES),1),...
            string({WFILES.name})');

disp(FILES)

%==========================================================================
%% LOAD THE CHROMOSOME TABLES FROM THEIR MAT FILES
%==========================================================================
clearvars -except P FILES


MATDAT = cell(numel(FILES),1);

for i = 1:numel(FILES)

    MATDAT{i} = load(FILES(i));

end





%==========================================================================
%% CONCATINATE THE CHR TABLES
%==========================================================================
clearvars -except P FILES MATDAT




GNOM = MATDAT{1}.TAB;

for i = 2:numel(MATDAT)

    GNOM = [GNOM; MATDAT{i}.TAB];

end








%==========================================================================
%% EXPORT TABLE
%==========================================================================
clearvars -except P FILES MATDAT GNOM


save([P.mat P.f 'PGS_STEP020_OUTPUT.mat'],'GNOM');









