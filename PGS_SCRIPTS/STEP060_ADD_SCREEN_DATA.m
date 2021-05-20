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


PGS = load([P.mat P.f 'PGS_STEP058_OUTPUT.mat']);


GNOMAD  = PGS.GNOMAD;



Nreport(GNOMAD)
%==========================================================================
%% LOAD LIST OF GENES INCLUDED IN CURRENT GENE SCREENING PANELS
%==========================================================================
clc; clearvars -except P PGS GNOMAD




P.datapath = [P.csv P.f 'panels' P.f 'ALL_PANELS.xlsx'];
P.opts     = detectImportOptions(P.datapath);
SCREENS    = readtable(P.datapath,P.opts);
SCREENS    = convertvars(SCREENS,@iscell,'string');





Nreport(GNOMAD)
%==========================================================================
%% MARK WHETHER GNOMAD GENES ARE IN SCREENING PANEL
%==========================================================================
clc; clearvars -except P PGS GNOMAD SCREENS




GNO = GNOMAD;

GNO.PANEL     = strings(height(GNO),1);
GNO.PANEL_REF = strings(height(GNO),1);
GNO.PANEL_DIS = strings(height(GNO),1);
GNO.PANEL_TF  = zeros(height(GNO),1);



[ai,bi] = ismember(GNO.GENE,SCREENS.GENE);


GNO.PANEL(ai)     = SCREENS.PANEL(bi(ai),:);
GNO.PANEL_REF(ai) = SCREENS.TRANSCRIPT(bi(ai),:);
GNO.PANEL_DIS(ai) = SCREENS.DISORDER(bi(ai),:);
GNO.PANEL_TF      = ~(GNO.PANEL == "");






UG = unique(GNO.GENE(GNO.PANEL == "INVITAE CORE")); NG = numel(UG);
disp(' '); disp('INVITAE CORE'); disp(UG); fprintf('\nN: %0.f\n\n',NG);

UG = unique(GNO.GENE(GNO.PANEL == "INVITAE MAIN")); NG = numel(UG);
disp(' '); disp('INVITAE MAIN'); disp(UG); fprintf('\nN: %0.f\n\n',NG);

UG = unique(GNO.GENE(GNO.PANEL == "INVITAE EXTENDED")); NG = numel(UG);
disp(' '); disp('INVITAE EXTENDED'); disp(UG); fprintf('\nN: %0.f\n\n',NG);

UG = unique(GNO.GENE(GNO.PANEL == "MYRIAD")); NG = numel(UG);
disp(' '); disp('MYRIAD UNIQUE'); disp(UG); fprintf('\nN: %0.f\n\n',NG);


GNOMAD = GNO;


Nreport(GNOMAD)






%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD SCREENS


DATAS = PGS.DATAS;

DATAS.SCREENS = SCREENS;

DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP060_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');


