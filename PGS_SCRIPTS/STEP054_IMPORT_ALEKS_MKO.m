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


PGS = load([P.mat P.f 'PGS_STEP053_2_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;   % GNOMAD DATASET
DAW     = PGS.DAW;      % DAWES DATASETS
GUO     = PGS.GUO;      % GUO & GREGG (2019)
MKO     = PGS.MKO;      % MOUSE KO DATASET





%==========================================================================
%% IMPORT MOUSE KNOCKOUT DATA
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO





MKO1 = readtable([P.csv P.f 'aleks' P.f 'ALEKS_MKO1.xlsx'], 'Sheet','MAIN');
MKO2 = readtable([P.csv P.f 'aleks' P.f 'ALEKS_MKO2.xlsx'], 'Sheet','MAIN');
MKO3 = readtable([P.csv P.f 'aleks' P.f 'ALEKS_MKO3.xlsx'], 'Sheet','MAIN');



MKO1 = convertvars(MKO1,@iscell,'string');
MKO2 = convertvars(MKO2,@iscell,'string');
MKO3 = convertvars(MKO3,@iscell,'string');





%==========================================================================
%% MK04 COMBINE DATA
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO MKO1 MKO2 MKO3


MKO4A = readtable([P.csv P.f 'aleks' P.f 'ALEKS_MKO4.xlsx'], 'Sheet','MAIN');
MKO4B = readtable([P.csv P.f 'aleks' P.f 'ALEKS_MKO4.xlsx'], 'Sheet','MAIN_MKO');

MKO4A = convertvars(MKO4A,@iscell,'string');
MKO4B = convertvars(MKO4B,@iscell,'string');


[C,ia,ic] = unique(MKO4B.MOUSE_GENE,'stable');
MKO4B = MKO4B(ia,:);

i = sum(MKO4A.MM_GENE == MKO4B.MOUSE_GENE' ,2)>0;
MKO4A = MKO4A(i,:);

i = sum(MKO4B.MOUSE_GENE == MKO4A.MM_GENE' ,2)>0;
MKO4B = MKO4B(i,:);

MKO4A = sortrows(MKO4A,'MM_GENE');
MKO4B = sortrows(MKO4B,'MOUSE_GENE');

MKO4 = [MKO4A MKO4B];


%==========================================================================
%% MK01 CONVERT DATA TYPES
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO MKO1 MKO2 MKO3 MKO4


MK = MKO1;

peek(MK)

%MK.Properties.VariableNames{'area3_main_phenotype'} = 'mouselethal';
%MK.Properties.VariableNames{'MIMNumber'} = 'mim';
MK = movevars(MK,{'mim','mouselethal'},'After','mko_gene');
unique(MK.mouselethal)

MK.mouselethal(MK.mouselethal == "NA") = "NaN";
MK.mouselethal(MK.mouselethal == "nonlethal") = "0";
MK.mouselethal(MK.mouselethal == "lethal") = "1";
unique(MK.mouselethal)

TF = ismissing(MK,"NA");
MK = fillmissing(MK,'constant',"0",'DataVariables',@isstring,'MissingLocations',TF);

TF = ismissing(MK,"yes");
MK = fillmissing(MK,'constant',"1",'DataVariables',@isstring,'MissingLocations',TF);

MK.cell_essential(MK.cell_essential == "essential") = "1";


MK.humanlethal = double(MK.humanlethal);
MK.mouselethal = double(MK.mouselethal);
MK.cell_essential = double(MK.cell_essential);
MK.biallelic_lof_tolerant = double(MK.biallelic_lof_tolerant);
MK.area1_high_embyronic_expression = double(MK.area1_high_embyronic_expression);
MK.mim = double(MK.mim);

MK.mousetarget = (MK.mim == 0 & MK.mouselethal == 1);
MK = movevars(MK,{'mousetarget'},'Before','mko_gene');

MKO1 = MK;





%==========================================================================
%% MK02 CONVERT DATA TYPES
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO MKO1 MKO2 MKO3 MKO4


MK = MKO2;

peek(MK)

%MK.Properties.VariableNames{'area3_main_phenotype'} = 'mouselethal';
%MK.Properties.VariableNames{'MIMNumber'} = 'mim';
MK = movevars(MK,{'mim','mouselethal'},'After','mko_gene');
unique(MK.mouselethal)

MK.mouselethal(MK.mouselethal == "NA") = "NaN";
MK.mouselethal(MK.mouselethal == "nonlethal") = "0";
MK.mouselethal(MK.mouselethal == "lethal") = "1";
unique(MK.mouselethal)

TF = ismissing(MK,"NA");
MK = fillmissing(MK,'constant',"0",'DataVariables',@isstring,'MissingLocations',TF);

TF = ismissing(MK,"yes");
MK = fillmissing(MK,'constant',"1",'DataVariables',@isstring,'MissingLocations',TF);

MK.cell_essential(MK.cell_essential == "essential") = "1";


MK.humanlethal = double(MK.humanlethal);
MK.mouselethal = double(MK.mouselethal);
MK.cell_essential = double(MK.cell_essential);
MK.biallelic_lof_tolerant = double(MK.biallelic_lof_tolerant);
MK.area1_high_embyronic_expression = double(MK.area1_high_embyronic_expression);
MK.mim = double(MK.mim);

MK.mousetarget = (MK.mim == 0 & MK.mouselethal == 1);
MK = movevars(MK,{'mousetarget'},'Before','mko_gene');

MKO2 = MK;








%==========================================================================
%% ANNOTATE GNOMAD DATASET WITH MKO INFO
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO MKO1 MKO2 MKO3 MKO4


peek(GNOMAD)




MK = MKO1;
GNO = GNOMAD;

idx = sum(GNO.GENE == MK.mko_gene' ,2)>0;
[ai,bi] = ismember(GNO.GENE,MK.mko_gene);

VNAMES = MK.Properties.VariableNames;



h = height(GNO);

for i = 1:numel(VNAMES)
    if isa(MK.(VNAMES{i}),'string')

        GNO.(VNAMES{i}) = strings(h,1);
        GNO( ai , VNAMES{i} ) = MK(bi(ai),VNAMES{i});

    else

        GNO.(VNAMES{i}) = nan(h,1);
        GNO( ai , VNAMES{i} ) = MK(bi(ai),VNAMES{i});

    end
end







GNOMAD = GNO;

peek(GNOMAD)



% unique(GNOMAD.DGENE)

% GNOMAD COLUMNS   1-112
% DAWES  COLUMNS 113-175
% MKO    COLUMNS 176-208
%==========================================================================
%% ANNOTATE GNOMAD DATASET WITH MKO INFO
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO MKO1 MKO2 MKO3 MKO4


peek(GNOMAD)




MK = MKO4;
GNO = GNOMAD;

idx = sum(GNO.GENE == MK.HS_GENE' ,2)>0;
[ai,bi] = ismember(GNO.GENE,MK.HS_GENE);

VNAMES = MK.Properties.VariableNames;



h = height(GNO);

for i = 1:numel(VNAMES)
    if isa(MK.(VNAMES{i}),'string')

        GNO.(VNAMES{i}) = strings(h,1);
        GNO( ai , VNAMES{i} ) = MK(bi(ai),VNAMES{i});

    else

        GNO.(VNAMES{i}) = nan(h,1);
        GNO( ai , VNAMES{i} ) = MK(bi(ai),VNAMES{i});

    end
end



GNOMAD = GNO;

peek(GNOMAD)



% unique(GNOMAD.HS_GENE)
% unique(GNOMAD.MOUSE_CHR(~isnan(GNOMAD.MOUSE_CHR)))

% GNOMAD COLUMNS   1-112
% DAWES  COLUMNS 113-175
% MKO1   COLUMNS 176-208
% MKO2   COLUMNS 209-269
%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P GNOMAD DAW GUO MKO MKO1 MKO2 MKO3 MKO4


MKO = MKO1;
MK2 = MKO4;


save([P.mat P.f 'PGS_STEP054_OUTPUT.mat'],'GNOMAD','DAW','MKO','GUO','MK2');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP054_OUTPUT.csv'])



