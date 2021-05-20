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




PGS = load([P.mat P.f 'PGS_STEP053_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;

DAW     = PGS.DAWES;









%==========================================================================
%% ANNOTATE GNOMAD DATASET WITH DAWES INFO
%==========================================================================
clc; clearvars -except P PGS GNOMAD DAW


DAW.DAWES2_S1.Properties.VariableNames{'gene'} = 'DGENE';


DAWES = DAW.DAWES2_S1;

idx = sum(GNOMAD.GENE == DAWES.DGENE' ,2)>0;
[ai,bi] = ismember(GNOMAD.GENE,DAWES.DGENE);

VNAMES = DAWES.Properties.VariableNames;

GNO = GNOMAD;

h = height(GNO);

for i = 1:numel(VNAMES)
    if isa(DAWES.(VNAMES{i}),'string')

        GNO.(VNAMES{i}) = strings(h,1);
        GNO( ai , VNAMES{i} ) = DAWES(bi(ai),VNAMES{i});

    else

        GNO.(VNAMES{i}) = nan(h,1);
        GNO( ai , VNAMES{i} ) = DAWES(bi(ai),VNAMES{i});

    end
end


GNO.isdawes = zeros(height(GNO),1);
GNO.isdawes(  GNO.DGENE ~= ""  ) = 1;


GNOMAD = GNO;







%==========================================================================
%% FORMAT TABLE DATA
%==========================================================================
clc; clearvars -except P PGS GNOMAD DAW


peek(GNOMAD)



GNO = GNOMAD;
TF = ismissing(GNO,"Y");
GNO = fillmissing(GNO,'constant',"1",'DataVariables',@isstring,'MissingLocations',TF);
TF = ismissing(GNO,"N");
GNO = fillmissing(GNO,'constant',"0",'DataVariables',@isstring,'MissingLocations',TF);



GNO.omim                = double(GNO.omim);
GNO.mim_number          = double(GNO.mim_number);
GNO.human_lethal_A      = double(GNO.human_lethal_A);
GNO.human_lethal_B      = double(GNO.human_lethal_B);
GNO.constrained         = double(GNO.constrained);
GNO.lethal_MGI          = double(GNO.lethal_MGI);
GNO.mouse_ko            = double(GNO.mouse_ko);
GNO.lethal_mouse        = double(GNO.lethal_mouse);
GNO.cell_ko             = double(GNO.cell_ko);
GNO.any_constraint      = double(GNO.any_constraint);
GNO.any_reg_constraint  = double(GNO.any_reg_constraint);
GNO.lethal_IMPC         = double(GNO.lethal_IMPC);
GNO.IMPC_ko             = double(GNO.IMPC_ko);
GNO.lethal_het_MGI      = double(GNO.lethal_het_MGI);
GNO.lethal_het_IMPC     = double(GNO.lethal_het_IMPC);
GNO.IMPC_het_ko         = double(GNO.IMPC_het_ko);
GNO.mouse_het_ko        = double(GNO.mouse_het_ko);
GNO.lethal_het_mouse    = double(GNO.lethal_het_mouse);
GNO.cell_essential      = double(GNO.cell_essential);


peek(GNO)

GNOMAD = GNO;





%==========================================================================
%% RENAME DAWES VARIABLES
%==========================================================================
clc; clearvars -except P PGS GNOMAD DAW



GNO = GNOMAD; peek(GNO)


GNO = movevars(GNO,{'DGENE','mim_number','omim'},'Before','hgnc_id');

GNO.Properties.VariableNames('DGENE') = {'GENES'};


unique(GNO.mim_number(~isnan(GNO.mim_number)))
unique(GNO.omim(~isnan(GNO.omim)))



GNO.Properties.VariableNames('omim') = {'has_mim'};
GNO.Properties.VariableNames('mim_number') = {'mim'};

GNO.has_mim(~isnan(GNO.mim)) = 1;
GNO.has_mim(isnan(GNO.mim))  = 0;

unique(GNO.has_mim)




VNAMES = string(GNO.Properties.VariableNames');

DAW_START = find(VNAMES == "GENES");
DAW_END   = find(VNAMES == "isdawes");


for i = 1:(DAW_END-DAW_START+1)

    j = DAW_START+i-1;

    VNAMES(j) = strcat("DAW_",VNAMES(j));

end

GNO.Properties.VariableNames = VNAMES;



GNOMAD = GNO;





%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD DAW


DAWES = DAW;

DATAS = PGS.DATAS;

save([P.mat P.f 'PGS_STEP056_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');


