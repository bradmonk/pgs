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


PGS = load([P.mat P.f 'PGS_STEP060_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;





%==========================================================================
%% IMPORT pLI INFO
%==========================================================================
clc; clearvars -except P PGS GNOMAD




P.lof_metrics = [P.csv P.f 'gnomad/gnomad_downloads/lof_metrics_by_gene.csv'];
OPT = detectImportOptions(P.lof_metrics);
LOFMET = readtable(P.lof_metrics,OPT);
LOFMET = convertvars(LOFMET,@iscell,'string');






%==========================================================================
%% ATTACH pLI INFO TO GNOMAD TABLE
%==========================================================================
clc; clearvars -except P PGS GNOMAD LOFMET





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


GNO = movevars(GNO,{'obs_lof','exp_lof','oe_lof'},'Before','DAW_GENES');
GNO = movevars(GNO,{'obs_mis','exp_mis','oe_mis'},'Before','DAW_GENES');
GNO = movevars(GNO,{'obs_syn','exp_syn','oe_syn'},'Before','DAW_GENES');
GNO = movevars(GNO,{'pLI'},'Before','DAW_GENES');


GNOMAD = GNO;

GNOMAD = GENEijk(GNOMAD);






%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD SCREENS LOFMET


DATAS = PGS.DATAS;

DATAS.LOFMET = LOFMET;

DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP070_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');





