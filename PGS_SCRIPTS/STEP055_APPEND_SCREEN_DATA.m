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
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP054_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;   % GNOMAD DATASET
DAW     = PGS.DAW;      % DAWES (2019)
GUO     = PGS.GUO;      % GUO & GREGG (2019)
MKO     = PGS.MKO;      % MOUSE KO DATASET SET-1
MK2     = PGS.MK2;      % MOUSE KO DATASET SET-2






Nreport(GNOMAD)
%==========================================================================
%% LOAD LIST OF GENES INCLUDED IN CURRENT GENE SCREENING PANELS
%==========================================================================
clc; clearvars -except P GNOMAD DAW MKO MK2 GUO




P.datapath = [P.csv P.f 'panels' P.f 'invitae_panel_main.xlsx'];
P.opts     = detectImportOptions(P.datapath);
INVITAE    = readtable(P.datapath,P.opts);
INVITAE.GENES = string(INVITAE.GENES);



P.datapath = [P.csv P.f 'panels' P.f 'invitae_panel_extended.xlsx'];
P.opts     = detectImportOptions(P.datapath);
INVITAX    = readtable(P.datapath,P.opts);
INVITAX.GENES = string(INVITAX.GENES);


P.datapath = [P.csv P.f 'panels' P.f 'myriad_panel_main.xlsx'];
P.opts     = detectImportOptions(P.datapath);
MYRIAD     = readtable(P.datapath,P.opts);
MYRIAD.GENES = string(MYRIAD.GENES);




SCREENS.INVITAE = INVITAE;
SCREENS.INVITAX = INVITAX;
SCREENS.MYRIAD  = MYRIAD;





Nreport(INVITAE, INVITAX, MYRIAD);
%==========================================================================
%% MARK WHETHER GNOMAD GENES ARE IN SCREENING PANEL
%==========================================================================
clc; clearvars -except P GNOMAD DAW MKO MK2 GUO SCREENS




% MARK INVITAE PANEL GENES
GNO = GNOMAD;
Gi = sum(GNO.GENE == SCREENS.INVITAE.GENES' ,2)>0;
GNO.INVITAE = zeros(height(GNO),1);
GNO.INVITAE(Gi) = 1;
GNOMAD = GNO;


% MARK INVITAX PANEL GENES
GNO = GNOMAD;
Gi = sum(GNO.GENE == SCREENS.INVITAX.GENES' ,2)>0;
GNO.INVITAX = zeros(height(GNO),1);
GNO.INVITAX(Gi) = 1;
GNOMAD = GNO;


% MARK MYRIAD PANEL GENES
GNO = GNOMAD;
Gi = sum(GNO.GENE == SCREENS.MYRIAD.GENES' ,2)>0;
GNO.MYRIAD = zeros(height(GNO),1);
GNO.MYRIAD(Gi) = 1;
GNOMAD = GNO;





numel(unique(GNOMAD.GENE(GNOMAD.INVITAE == 1)))
numel(unique(GNOMAD.GENE(GNOMAD.INVITAX == 1)))
numel(unique(GNOMAD.GENE(GNOMAD.MYRIAD == 1)))



Nreport(GNOMAD)
%==========================================================================
%% IMPORT pLI INFO
%==========================================================================
clc; clearvars -except P GNOMAD DAW MKO MK2 GUO




P.lof_metrics = [P.csv P.f 'gnomad/gnomad_downloads/lof_metrics_by_gene.csv'];
OPT = detectImportOptions(P.lof_metrics);
LOFMET = readtable(P.lof_metrics,OPT);
LOFMET = convertvars(LOFMET,@iscell,'string');






%==========================================================================
%% ATTACH pLI INFO TO GNOMAD TABLE
%==========================================================================
clc; clearvars -except P GNOMAD DAW MKO MK2 GUO LOFMET





GNOMAD = GENEijk(GNOMAD);





GNO = GNOMAD;

idx = sum(GNOMAD.GENE == LOFMET.gene' ,2)>0;
[ai,bi] = ismember(GNO.GENE,LOFMET.gene);

GNO.obs_syn(ai) = LOFMET.obs_syn(bi(ai));
GNO.exp_syn(ai) = LOFMET.exp_syn(bi(ai));
GNO.oe_syn(ai)  = LOFMET.oe_syn(bi(ai));
GNO.obs_mis(ai) = LOFMET.obs_mis(bi(ai));
GNO.exp_mis(ai) = LOFMET.exp_mis(bi(ai));
GNO.oe_mis(ai)  = LOFMET.oe_mis(bi(ai));
GNO.obs_lof(ai) = LOFMET.obs_lof(bi(ai));
GNO.exp_lof(ai) = LOFMET.exp_lof(bi(ai));
GNO.oe_lof(ai)  = LOFMET.oe_lof(bi(ai));
GNO.pLI(ai)     = LOFMET.pLI(bi(ai));


GNO = movevars(GNO,{'obs_lof','exp_lof','oe_lof'},'After','isSV');
GNO = movevars(GNO,{'obs_mis','exp_mis','oe_mis'},'After','isSV');
GNO = movevars(GNO,{'obs_syn','exp_syn','oe_syn'},'After','isSV');
GNO = movevars(GNO,{'pLI'},'After','isSV');


GNOMAD = GNO;

GNOMAD = GENEijk(GNOMAD);





% GNOMAD COLUMNS   1-114
% PLI    COLUMNS 115-124
% DAWES  COLUMNS 125-186
% MKO1   COLUMNS 187-219
% MKO2   COLUMNS 220-280
% SCR    COLUMNS 281-283
%==========================================================================
%% CREATE MK2 TARGET COLUMN
%==========================================================================
clc; clearvars -except P GNOMAD DAW MKO MK2 GUO LOFMET




GNOMAD.Properties.VariableNames{'mousetarget'} = 'MKO1_TARGETS';
GNOMAD.Properties.VariableNames{'VIABILITY_CONFIRMED'} = 'MKO2_TARGETS';
GNOMAD.Properties.VariableNames{'CANDIDATE_GENE'} = 'DAWES_TARGETS';



GNOMAD.MKO2_TARGETS = (GNOMAD.VIABILITY == "Lethal") & (GNOMAD.omim ~= "Y");

MK1 = MKO;


%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P GNOMAD DAW MK1 MK2 GUO


save([P.mat P.f 'PGS_STEP055_OUTPUT.mat'],'GNOMAD','DAW','MK1','MK2','GUO');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP054_OUTPUT.csv'])



