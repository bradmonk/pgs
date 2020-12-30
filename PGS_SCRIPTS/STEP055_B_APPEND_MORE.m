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


PGS = load([P.mat P.f 'PGS_STEP055_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;   % GNOMAD DATASET
DAW     = PGS.DAW;      % DAWES (2019)
GUO     = PGS.GUO;      % GUO & GREGG (2019)
MK1     = PGS.MK1;      % MOUSE KO DATASET SET-1
MK2     = PGS.MK2;      % MOUSE KO DATASET SET-2








peek(GNOMAD)


% GNOMAD COLUMNS   1-114
% PLI    COLUMNS 115-124
% DAWES  COLUMNS 125-186
% MKO1   COLUMNS 187-219
% MKO2   COLUMNS 220-280
% SCR    COLUMNS 281-283
%==========================================================================
%% IMPORT MOUSE KNOCKOUT DATA
%==========================================================================
clc; clearvars -except P GNOMAD DAW MK1 MK2 GUO





TBL1 = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','HUMAN_LETHAL');
TBL2 = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','MOUSE_LETHAL_A');
TBL3 = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','MOUSE_LETHAL_B');
TBL4 = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','HOMOLOGUES');



TBL1 = convertvars(TBL1,@iscell,'string');
TBL2 = convertvars(TBL2,@iscell,'string');
TBL3 = convertvars(TBL3,@iscell,'string');
TBL4 = convertvars(TBL4,@iscell,'string');






%==========================================================================
%% FORMAT TABLE DATA
%==========================================================================
clc; clearvars -except P GNOMAD DAW MK1 MK2 GUO TBL1 TBL2 TBL3 TBL4
%----------------------------------------
% GNOMAD COLUMNS   1-114
% PLI    COLUMNS 115-124
% DAWES  COLUMNS 125-186
% MKO1   COLUMNS 187-219
% MKO2   COLUMNS 220-280
% SCR    COLUMNS 281-283
% 
% A GOOD PLACE TO LOOK AT DAWES DATA IS AT MUTYH & TOE1 GENES
% 
%   GTAB_MUTYH = GNOMAD(GNOMAD.GENE == "MUTYH",:);
%   GTAB_TOE1  = GNOMAD(GNOMAD.GENE == "TOE1",:);
% 
%----------------------------------------
peek(GNOMAD)



GNO = GNOMAD;
TF = ismissing(GNO,"Y");
GNO = fillmissing(GNO,'constant',"1",'DataVariables',@isstring,'MissingLocations',TF);
TF = ismissing(GNO,"N");
GNO = fillmissing(GNO,'constant',"0",'DataVariables',@isstring,'MissingLocations',TF);




GNO.omim            = double(GNO.omim);
GNO.mim_number      = double(GNO.mim_number);
GNO.human_lethal_A  = double(GNO.human_lethal_A);
GNO.human_lethal_B  = double(GNO.human_lethal_B);
GNO.constrained     = double(GNO.constrained);
GNO.lethal_MGI      = double(GNO.lethal_MGI);
GNO.mouse_ko        = double(GNO.mouse_ko);
GNO.lethal_mouse    = double(GNO.lethal_mouse);
GNO.cell_ko         = double(GNO.cell_ko);
%GNO.cell_essential_hits
%GNO.cell_essential
%GNO.DAWES_TARGETS
%GNO.isdawes



% OMIM COLUMNS
% GNO.omim
% GNO.mim_number
% GNO.mim

GNOMAD = GNO;




%return
%==========================================================================
%% FLAG TARGET GENES
%==========================================================================
clc; clearvars -except P GNOMAD DAW MK1 MK2 GUO TBL1 TBL2 TBL3 TBL4



peek(GNOMAD)



GNO = GNOMAD;




i = sum(GNO.GENE == TBL1.HL_gene' ,2)>0;
GNO.HUMAN_TARGETS = zeros(height(GNO),1);
GNO.HUMAN_TARGETS(i) = 1;




i = sum(GNO.GENE == TBL2.ML1_gene' ,2)>0;
GNO.MOUSE1_GENES      = zeros(height(GNO),1);
GNO.MOUSE1_TARGETS    = zeros(height(GNO),1);
GNO.MOUSE1_GENES(i)   = 1;
GNO.MOUSE1_TARGETS(i) = 1;



i = sum(GNO.GENE == TBL3.ML2_hs_gene' ,2)>0;
GNO.MOUSE2_GENES      = zeros(height(GNO),1);
GNO.MOUSE2_TARGETS    = zeros(height(GNO),1);
GNO.MOUSE2_GENES(i)   = 1;
GNO.MOUSE2_TARGETS(i) = 1;





% UPDATED MOUSE2_TARGETS RAW DATA IN EXCEL SHEET USING...
%----------------------------------------
% i = sum(TBL3.ML2_gene == TBL4.mm_gene' ,2)>0;
% j = sum(TBL4.mm_gene == TBL3.ML2_gene' ,2)>0;
% T3 = TBL3(i,:);
% T4 = TBL4(j,:);
% TBL3.ML2_homologous = i;
% TBL3.ML2_homologous = i;
%----------------------------------------
% 
% 
% BUT....
% AFTER FINDING MATCHES IN GNOMAD TABLE, THIS SET
% SHOULD BE FILTERED TO EXCLUDE GENES WITH MIM NUMBERS
% 
% OMIM COLUMNS IN GNOMAD TABLE:
%   GNOMAD.omim         (COL 126)
%   GNOMAD.mim_number   (COL 127)
%   GNOMAD.mim          (COL 186)
%----------------------------------------
% 
% 
% 
% BUT...
% FIRST FIX DISCREPANCIES BETWEEN GNOMAD OMIM DATA
%-----------------------
% MIM.HAS_OMIM(151) IS MARKED "0"
% MIM.OMIM1(151) DOESNT HAVE MIM NUMBER
% MIM.OMIM2(151) DOES   HAVE MIM NUMBER
%----------------------------------------


HAS_OMIM = GNO.omim;
OMIM1    = GNO.mim_number;
OMIM2    = GNO.mim;

MIM = table(HAS_OMIM,OMIM1,OMIM2);

MIM.HAS_OMIM(isnan(MIM.HAS_OMIM)) = 0;
MIM.OMIM1(isnan(MIM.OMIM1)) = 0;
MIM.OMIM2(isnan(MIM.OMIM2)) = 0;

MIM.HAS_OMIM = (MIM.OMIM1 > 0) | (MIM.OMIM2 > 0);

GNO = [GNO , MIM];

GNO.MOUSE2_TARGETS(GNO.HAS_OMIM==1) = 0;







[C,ia,ic] = unique(GNO.GENE,'stable');
TAB = GNO(ia,:);
TAB.DAWES_TARGETS = TAB.DAWES_TARGETS == 1;
clc
fprintf('GNOMAD LoF GENES:	%9.0f \n' , height(TAB))
fprintf('~DAWES:            %10.0f \n' , sum(~TAB.DAWES_TARGETS))
fprintf(' DAWES:            %10.0f \n' , sum(TAB.DAWES_TARGETS))
fprintf(' HUMAN:            %10.0f \n' , sum(TAB.HUMAN_TARGETS))
fprintf(' MOUSE1:           %10.0f \n' , sum(TAB.MOUSE1_TARGETS))
fprintf(' MOUSE2:           %10.0f \n' , sum(TAB.MOUSE2_TARGETS))
fprintf(' DAWES & HUMAN:    %10.0f \n' , sum(TAB.DAWES_TARGETS & TAB.HUMAN_TARGETS))
fprintf(' DAWES & MOUSE1:   %10.0f \n' , sum(TAB.DAWES_TARGETS & TAB.MOUSE1_TARGETS))
fprintf(' DAWES & MOUSE2:   %10.0f \n' , sum(TAB.DAWES_TARGETS & TAB.MOUSE2_TARGETS))
fprintf(' HUMAN & MOUSE1:   %10.0f \n' , sum(TAB.HUMAN_TARGETS & TAB.MOUSE1_TARGETS))
fprintf(' HUMAN & MOUSE2:   %10.0f \n' , sum(TAB.HUMAN_TARGETS & TAB.MOUSE2_TARGETS))
sum(TAB.DAWES_TARGETS & TAB.HAS_OMIM)
sum(TAB.HUMAN_TARGETS & TAB.HAS_OMIM)
sum(TAB.MOUSE1_TARGETS & TAB.HAS_OMIM)
sum(TAB.MOUSE2_TARGETS & TAB.HAS_OMIM)





GNO = movevars(GNO,{'DAWES_TARGETS','HUMAN_TARGETS','MOUSE1_TARGETS','MOUSE2_TARGETS'},'After','OMIM2');
GNO.MOUSE1_GENES = [];
GNO.MOUSE2_GENES = [];



GNOMAD = GNO;

peek(GNOMAD)





%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P GNOMAD DAW MK1 MK2 GUO TBL1 TBL2 TBL3 TBL4


TABLES.DAW = DAW;
TABLES.MK1 = MK1;
TABLES.MK2 = MK2;
TABLES.GUO = GUO;
TABLES.HUM = TBL1;
TABLES.MS1 = TBL2;
TABLES.MS2 = TBL3;
TABLES.OLO = TBL4;

save([P.mat P.f 'PGS_STEP055_B_OUTPUT.mat'],'GNOMAD','TABLES');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP055_B_OUTPUT.csv'])



