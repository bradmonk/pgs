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




PGS = load([P.mat P.f 'PGS_STEP053_1_CONCAT_STRUCTVARS.mat']);
%PGS2 = load([P.mat P.f 'PGS_STEP052_OUTPUT.mat']);



GNOMAD  = PGS.GNOMAD;	% GNOMAD DATASET
GUO     = PGS.GUO;      % GUO & GREGG (2019)
MKO     = PGS.MKO;      % MOUSE KO DATASET
DAW     = PGS.DAW;      % DAWES_S2_S3.mat


%DAW = load([P.mat P.f 'DAWES_S2_S3.mat']);






%==========================================================================
%% FLAG THE 3435 CANDIDATE DEVELOPMENTAL LETHAL HUMAN GENES
%==========================================================================
clc; clearvars -except P PGS GNOMAD GUO MKO DAW



% CANDIDATE DEVELOPMENTAL LETHAL HUMAN GENES (3435 TOTAL)
%--------------------------------------------------------------------------
% We therefore curate 3435 ‘candidate developmental lethal’ human genes:
% essential for murine development or cellular viability, not yet linked to
% human disorders, presenting strong candidates for unexplained infertility
% and prenatal/infantile mortality.
% 
% 
% GENETIC CONSTRAINT
% Adaptations will often be imperfect because of genetic constraints. An
% example of such a constraint occurs when the heterozygote at a locus has
% a higher fitness than either homozygote, and the population evolves to an
% equilibrium at which all three genotypes are present. A proportion of the
% population must therefore have the deleterious homozygous genotypes.
% 
% This situation arises because the heterozygotes cannot, under Mendelian
% inheritance, produce purely heterozygous offspring: they cannot breed
% true. Insofar as heterozygous advantage exists, some members of the
% natural populations will be imperfectly adapted. The importance of
% heterozygous advantage is controversial; but there are undoubted
% examples, such as sickle cell anemia - caused by sickle shaped blood
% cells, pictured opposite - which is a practical example of imperfect
% adaptation due to a genetic constraint.
%--------------------------------------------------------------------------





i = sum(DAW.DAWES2_S1.DGENE == DAW.DAWES2_S3.GENES' ,2)>0;
DAW.DAWES2_S1.CANDIDATE_GENE = i;










%==========================================================================
%% CREATE GNOMAD DATA TABLE THAT INCLUDES ONLY DAWES GENES
%==========================================================================
%{

clc; clearvars -except P PGS GNOMAD GUO MKO DAW



DAWES = DAW.DAWES2_S1;

idx = sum(GNOMAD.GENE == DAWES.GENE ,2)>0;

DAWMAD = GNOMAD(idx , :);






h = height(DAWMAD);
DAWMAD.mim         = zeros(h,1);
DAWMAD.geneset     = strings(h,1);
DAWMAD.autosomal   = zeros(h,1) == 1;
DAWMAD.somatic     = zeros(h,1) == 1;
DAWMAD.xlinked     = zeros(h,1) == 1;
DAWMAD.recessive   = zeros(h,1) == 1;
DAWMAD.mosaic      = zeros(h,1) == 1;
DAWMAD.digenic     = zeros(h,1) == 1;
DAWMAD.mito        = zeros(h,1) == 1;
DAWMAD.inheritance = strings(h,1);
DAWMAD.matches     = strings(h,1);


[ai,bi] = ismember(DAWMAD.GENE,DAWES.GENE);



DAWMAD.mim(ai)          = DAWES.mim(bi(ai),:);
DAWMAD.geneset(ai)      = DAWES.gene(bi(ai),:);
DAWMAD.autosomal(ai)    = DAWES.autosomal(bi(ai),:);
DAWMAD.somatic(ai)      = DAWES.somatic(bi(ai),:);
DAWMAD.xlinked(ai)      = DAWES.xlinked(bi(ai),:);
DAWMAD.recessive(ai)    = DAWES.recessive(bi(ai),:);
DAWMAD.mosaic(ai)       = DAWES.mosaic(bi(ai),:);
DAWMAD.digenic(ai)      = DAWES.digenic(bi(ai),:);
DAWMAD.mito(ai)         = DAWES.mito(bi(ai),:);
DAWMAD.inheritance(ai)  = DAWES.inheritance(bi(ai),:);
DAWMAD.matches(ai)      = DAWES.matches(bi(ai),:);
%}




%==========================================================================
%% ANNOTATE GNOMAD DATASET WITH DAWES INFO
%==========================================================================
clc; clearvars -except P PGS GNOMAD GUO MKO DAW


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
%% ANNOTATE GNOMAD DATASET WITH DAWES INFO
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD GUO MKO DAW



h = height(GNOMAD);
GNOMAD.mim         = zeros(h,1);
GNOMAD.geneset     = strings(h,1);
GNOMAD.autosomal   = zeros(h,1) == 1;
GNOMAD.somatic     = zeros(h,1) == 1;
GNOMAD.xlinked     = zeros(h,1) == 1;
GNOMAD.recessive   = zeros(h,1) == 1;
GNOMAD.mosaic      = zeros(h,1) == 1;
GNOMAD.digenic     = zeros(h,1) == 1;
GNOMAD.mito        = zeros(h,1) == 1;
GNOMAD.inheritance = strings(h,1);
GNOMAD.matches     = strings(h,1);



[ai,bi] = ismember(GNOMAD.GENE,DAWMAD.GENE);


GNOMAD.mim(ai)          = DAWMAD.mim(bi(ai),:);
GNOMAD.geneset(ai)      = DAWMAD.geneset(bi(ai),:);
GNOMAD.autosomal(ai)    = DAWMAD.autosomal(bi(ai),:);
GNOMAD.somatic(ai)      = DAWMAD.somatic(bi(ai),:);
GNOMAD.xlinked(ai)      = DAWMAD.xlinked(bi(ai),:);
GNOMAD.recessive(ai)    = DAWMAD.recessive(bi(ai),:);
GNOMAD.mosaic(ai)       = DAWMAD.mosaic(bi(ai),:);
GNOMAD.digenic(ai)      = DAWMAD.digenic(bi(ai),:);
GNOMAD.mito(ai)         = DAWMAD.mito(bi(ai),:);
GNOMAD.inheritance(ai)  = DAWMAD.inheritance(bi(ai),:);
GNOMAD.matches(ai)      = DAWMAD.matches(bi(ai),:);



DAWMAD.isdawes   = zeros(height(DAWMAD),1) == 0;
GNOMAD.isdawes   = zeros(height(GNOMAD),1) == 1;
GNOMAD.isdawes(GNOMAD.mim > 0) = 1;
%}



%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD GUO MKO DAW


save([P.mat P.f 'PGS_STEP053_2_OUTPUT.mat'],'GNOMAD','MKO','GUO','DAW');





%==========================================================================
%% IMPORT OFFICIAL GENE NAMES FOR DAWES 2019 GENE LIST
%==========================================================================
%{
% clc; clearvars -except P PGS GNOMAD GUO MKO DAW
% 
% 
% 
% 
% HGNC = readtable([P.csv P.f 'HGNC_RESULTS_FOR_DAWES_GENES.csv']);
% %HGNC = readtable([P.csv P.f 'DAWES_HGNC_COMBINED_DATA.csv']);
% peek(HGNC)
% 
% 
% HGNC.Input          = string(HGNC.Input);
% HGNC.MatchType      = string(HGNC.MatchType);
% HGNC.ApprovedSymbol = string(HGNC.ApprovedSymbol);
% HGNC.ApprovedName   = string(HGNC.ApprovedName);
% HGNC.HGNCID         = string(HGNC.HGNCID);
% HGNC.Location       = string(HGNC.Location);
% 
% 
% HGNC.Properties.VariableNames{'ApprovedSymbol'} = 'GENE';
% 
% 
% % HGNC(HGNC.GENE == "",:) = [];
% % [C,ia,ic] = unique(HGNC.GENE,'stable');
% % HGNC = HGNC(ia,:);






%==========================================================================
%% MERGE DAWES INFO WITH HGNC INFO
%==========================================================================
% clc; clearvars -except P PGS GNOMAD GUO MKO DAW
% 
% GENES = DAWES.gene;
% 
% GENELIST = cellfun(@split, GENES,  repmat({";"}, size(GENES,1) ,1)  ,...
%     'UniformOutput',false );
% 
% nCells = cellfun(@numel, GENELIST);
% 
% rows = numel(nCells);
% cols = max(nCells);
% 
% GEN = strings(rows,cols);
% 
% for i = 1:rows
% 
%     s = string( GENELIST{i} )';
% 
%     GEN(i , 1:numel(s)) = s;
% 
% end
% 
% GLIST = GEN(:);
% GLIST = unique(GLIST);
% 
% if isempty(GLIST{1})
%     GLIST(1) = [];
% end
% GLIST(:)
% 
% 
% DAWES.GEN = GEN;
% 
% DAWES.GEN = strtrim(DAWES.GEN);
% DAWES.GEN = strip(DAWES.GEN);

% DAWGEN = GLIST;
% DAWGEN = HGNC.ApprovedSymbol(HGNC.ApprovedSymbol ~= "");




%==========================================================================
%% MERGE DAWES INFO WITH HGNC INFO
%==========================================================================
% clc; clearvars -except P PGS GNOMAD GUO MKO DAW
% 
% 
% 
% DAW = DAWES(:,[2 4:10]);
% 
% [C,ia,ic] = unique(DAW,'stable','rows');
% 
% DAW = DAWES(ia,:);
% DAW.gene_mim_nums = [];
% DAW.search_field = [];
% DAW.mim_entry_type = [];
% DAW.pheno_map_keys = [];
% 
% 
% 
% 
% 
% 
% DAW.gene        = string(DAW.gene);
% DAW.autosomal   = string(DAW.autosomal) == "TRUE";
% DAW.somatic     = string(DAW.somatic) == "TRUE";
% DAW.xlinked     = string(DAW.xlinked) == "TRUE";
% DAW.recessive   = string(DAW.recessive) == "TRUE";
% DAW.mosaic      = string(DAW.mosaic) == "TRUE";
% DAW.digenic     = string(DAW.digenic) == "TRUE";
% DAW.mito        = string(DAW.mito) == "TRUE";
% DAW.inheritance = string(DAW.inheritance);
% DAW.matches     = string(DAW.matches);
% 
% 
% 
% h = height(HGNC);
% HGNC.mim         = zeros(h,1);
% HGNC.gene        = strings(h,1);
% HGNC.autosomal   = zeros(h,1) == 1;
% HGNC.somatic     = zeros(h,1) == 1;
% HGNC.xlinked     = zeros(h,1) == 1;
% HGNC.recessive   = zeros(h,1) == 1;
% HGNC.mosaic      = zeros(h,1) == 1;
% HGNC.digenic     = zeros(h,1) == 1;
% HGNC.mito        = zeros(h,1) == 1;
% HGNC.inheritance = strings(h,1);
% HGNC.matches     = strings(h,1);
% 
% 
% 
% HGNC.GENE(HGNC.GENE == "") = "NA";
% 
% for i = 1:height(DAW)
% 
%     m = sum(HGNC.Input == DAW.GEN(i,:) , 2) > 0;
% 
%     HGNC.mim(m)         = DAW.mim(i);
%     HGNC.gene(m)        = DAW.gene(i);
%     HGNC.autosomal(m)   = DAW.autosomal(i);
%     HGNC.somatic(m)     = DAW.somatic(i);
%     HGNC.xlinked(m)     = DAW.xlinked(i);
%     HGNC.recessive(m)   = DAW.recessive(i);
%     HGNC.mosaic(m)      = DAW.mosaic(i);
%     HGNC.digenic(m)     = DAW.digenic(i);
%     HGNC.mito(m)        = DAW.mito(i);
%     HGNC.inheritance(m) = DAW.inheritance(i);
%     HGNC.matches(m)     = DAW.matches(i);
% 
% end
% 
% 
% DAWES = HGNC;
% 
% DAWES.Properties.VariableNames{'Input'} = 'DAWGENE';


%==========================================================================
%% IMPORT FULL HUMAN GENE LIST
%==========================================================================
% clc; clearvars -except P PGS GNOMAD GUO MKO DAW
% 
% 
% ALLG = readtable([P.csv P.f 'FULL_HUMAN_GENE_LIST.xlsx']);
% peek(ALLG)
% V = ALLG.Properties.VariableNames;
% for i = 1:numel(V)
%     ALLG.(V{i}) = string(ALLG.(V{i}));
% end
%}
