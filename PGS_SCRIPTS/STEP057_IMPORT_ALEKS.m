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


PGS = load([P.mat P.f 'PGS_STEP056_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;







%==========================================================================
%% IMPORT ALL GENES DATA SHEET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


PGS.DATAS.ALL = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','ALL_GENES');

PGS.DATAS.ALL = convertvars(PGS.DATAS.ALL,@iscell,'string');


%==========================================================================
%% HUMAN & MOUSE TARGETS
%==========================================================================
clc; clearvars -except P PGS GNOMAD





% TARGS.HUM = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','HUMAN_LETHAL');
% TARGS.MK1 = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','MOUSE_LETHAL_A');
% TARGS.MK2 = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','MOUSE_LETHAL_B');
% TARGS.HMS = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','HOMOLOGUES');
% TARGS.ALL = readtable([P.csv P.f 'aleks' P.f 'PGS_MAIN.xlsx'], 'Sheet','ALL_GENES');
% 
% 
% 
% TARGS.HUM = convertvars(TARGS.HUM,@iscell,'string');
% TARGS.MK1 = convertvars(TARGS.MK1,@iscell,'string');
% TARGS.MK2 = convertvars(TARGS.MK2,@iscell,'string');
% TARGS.HMS = convertvars(TARGS.HMS,@iscell,'string');
% TARGS.ALL = convertvars(TARGS.ALL,@iscell,'string');









%==========================================================================
%% INSTRUCTIONS
%==========================================================================
%{





STEP 1. BUIILD HUMAN LETHAL GENE LIST


    

    Use the human lethal genes from Dawes paper (n=624) containing genes 
    known to elicit lethality in humans (as described in OMIM).

    Filter the genes that are recessive (TARGS.HUM.HL_phenotypes; include
    "recessive" and "recessive/dominant" entries)



    TARGS.HUM                   Dawes raw dataset
    TARGS.HUM.HL_phenotypes     column of recessive/dominant phenotypes






STEP 2. BUILD MOUSE LETHAL GENE LIST-A


    Import the mouse lethal genes List A (placental genes master file)
    
    Keep only rows that contain "lethal" genes. The column labeled
    "ML1_area3_main_phenotype" contains which genes are lethal.

    Keep only genes that have a MIM Number. The column labeled 
    "ML1_mim_num" contains mim number entries. You should end 
    up with n=892 genes.



    TARGS.MK1                               mouse list-a raw data
    TARGS.MK1.ML1_area3_main_phenotype      column indicating "lethal" genes
    TARGS.MK1.ML1_mim_num                   column containing mim numbers





STEP 3. BUILD MOUSE LETHAL GENE LIST-B


    Import the mouse lethal genes List B (DMDD - early lethals)

    All the genes listed are lethal ("A line is classified as early
    lethal if no live homozygous embryos are recovered at E9.5").


    Determine human homologues. Column "ML2_hs_gene" contains the
    human homologe gene name, if one exists. Keep only the genes with 
    human homologues.

    Remove the genes without a MIM number. This table does not contain
    its own mim number column so this will need to be cross referenced 
    to a larger list of human genes. Use TARGS.ALL imported above as
    the cross reference between "MK2.ML2_hs_gene" and "ALL.gene"



    TARGS.MK2                               mouse list-b raw data
    TARGS.MK2.ML2_hs_gene                   human homologe gene names
    TARGS.ALL.gene                          column containing mim numbers




STEP 4. COMBINE MOUSE LIST-1 & LIST-2

    Merge mouse list-1 and list-2 into a single file

    Remove duplicates Take out the duplicates





%}















%==========================================================================
%% DETERMINE STARTING NUMBER OF HUMAN TARGET GENES
%==========================================================================
clc; clearvars -except P PGS GNOMAD TARGS


HUM = PGS.DATAS.HUMAN1;
%HUM = TARGS.HUM;


N.TOT = height(HUM);





HUM.HUM_OMIM        = HUM.HL_mim_num;
HUM.HUM_LETHAL      = HUM.HL_human_lethal == "yes";
HUM.IS_RECESSIVE    = HUM.HL_human_lethal == "yes";



% KEEP ONLY GENES WITH MIM NUMBERS (HUM.HL_mim_num)
%------------------------------------------
HUM.HL_OMIM = HUM.HL_mim_num;
HUM = HUM(HUM.HL_OMIM > 0,:);
N.MIM = height(HUM);





% KEEP ONLY LETHAL GENES (HUM.HL_human_lethal)
%------------------------------------------
HUM.HL_LETHAL = HUM.HL_human_lethal == "yes";
HUM = HUM(HUM.HL_LETHAL,:);
N.LTH = height(HUM);






% KEEP ONLY RECESSIVE GENES (HUM.HL_phenotypes)
%------------------------------------------
T = regexpi(HUM.HL_phenotypes,'recessive');
L(:,1) = ~(cellfun(@isempty,T));
HUM.HL_RECESSIVE = sum(L,2)>0;
HUM = HUM(HUM.HL_RECESSIVE == 1 ,:);
N.REC = height(HUM);






% REPORT SIZE OF TABLE AT EACH STEP
%------------------------------------------
N.TOT       % 624 TOTAL DAWES GENES
N.MIM       % 624 HAVE MIM NUMBERS
N.LTH       % 624 ALL ARE LETHAL
N.MIM       % 624 HAVE MIM NUMBERS
N.REC       % 527 RECESSIVE GENES



HUM = movevars(HUM,{'HL_OMIM','HL_LETHAL','HL_RECESSIVE'},'After','HL_gene');




%TARGS.HUM = HUM;
PGS.DATAS.HUMAN1 = HUM;



















%==========================================================================
%% DETERMINE STARTING NUMBER OF MOUSE TARGET LIST-A GENES
%==========================================================================
clc; clearvars -except P PGS GNOMAD TARGS



MK1 = PGS.DATAS.MOUSE1;
%MK1 = TARGS.MK1;




N.TOT = height(MK1);





% MK1.HAS_HOMOLOGUE
% MK1.IS_LETHAL
% MK1.NO_MIM_NUM





% CREATE A HOMOLOGES FLAG COLUMN
%------------------------------------------------------------
MK1.ML1_HAS_HOMOLOGUE = ones(height(MK1),1);
MK1 = MK1(MK1.ML1_HAS_HOMOLOGUE == 1 ,:);
N.HOG = height(MK1);






% ONLY USE GENES MARKED AS LETHAL
%------------------------------------------------------------
MK1.ML1_WHEN_DIED   = MK1.ML1_secondary_phenotype_cat;
MK1.ML1_LETHAL      = MK1.ML1_area3_main_phenotype;
MK1.ML1_ISLETHAL    = MK1.ML1_LETHAL == "lethal";
MK1 = MK1(MK1.ML1_ISLETHAL == 1 ,:);
N.LTH = height(MK1);






% KEEP ONLY GENES WITHOUT A MIM NUMBER
%------------------------------------------------------------
MK1.ML1_HASMIM = (MK1.ML1_mim_num ~= "NA");
MK1.ML1_NOMIM  = (MK1.ML1_mim_num == "NA");
MK1 = MK1(MK1.ML1_NOMIM , :);
N.MIM = height(MK1);






MK1 = movevars(MK1,{'ML1_HAS_HOMOLOGUE','ML1_ISLETHAL','ML1_LETHAL',...
    'ML1_WHEN_DIED','ML1_NOMIM','ML1_HASMIM'},'After','ML1_gene');


% ML1_HAS_HOMOLOGUE
% ML1_ISLETHAL
% ML1_LETHAL
% ML1_WHEN_DIED
% ML1_NOMIM
% ML1_HASMIM


N.TOT       % 892 TOTAL DAWES GENES
N.HOG       % 892 HOMOLOGES
N.LTH       % 892 ALL ARE LETHAL
N.MIM       % 892 HAVE MIM NUMBERS






% TARGS.MK1 = MK1;
PGS.DATAS.MOUSE1 = MK1;














%==========================================================================
%% DETERMINE STARTING NUMBER OF MOUSE TARGET LIST-B GENES
%==========================================================================
clc; clearvars -except P PGS GNOMAD TARGS



MK2 = PGS.DATAS.MOUSE2;
% MK2 = TARGS.MK2;


N.TOT = height(MK2);






% ONLY USE GENES MARKED AS LETHAL
%------------------------------------------------------------
MK2.ML2_ISLETHAL  = ones(height(MK2),1);
N.LTH = height(MK2);






% CREATE A HOMOLOGES FLAG COLUMN
%------------------------------------------------------------
MK2.ML2_HAS_HOMOLOGUE = (MK2.ML2_hs_gene ~= "");
MK2 = MK2(MK2.ML2_HAS_HOMOLOGUE==1,:);
N.HOG = height(MK2);







% KEEP ONLY GENES WITHOUT A MIM NUMBER
%------------------------------------------------------------
MK2.ML2_NOMIM  = ones(height(MK2),1);
MK2.ML2_HASMIM = zeros(height(MK2),1);
N.MIM = height(MK2);

% ALL = TARGS.ALL;
% i = (sum(ALL.gene == MK2.ML2_hs_gene' ,2)>0) & (ALL.omim ~= "Y");
% ALL = ALL(i,:);
% sum(MK2.ML2_hs_gene == ALL.gene' ,2)>0












MK2 = movevars(MK2,{'ML2_ISLETHAL','ML2_NOMIM','ML2_HASMIM','ML2_HAS_HOMOLOGUE'},'After','ML2_gene');



N.TOT       % 892 TOTAL DAWES GENES
N.LTH       % 892 ALL ARE LETHAL
N.MIM       % 892 HAVE MIM NUMBERS
N.HOG       % 892 HOMOLOGUES




% TARGS.MK2 = MK2;
PGS.DATAS.MOUSE2 = MK2;














% return
%==========================================================================
%% INSERT TARGET GENES INTO GNOMAD DATA
%==========================================================================
clc; clearvars -except P PGS GNOMAD


GNO = GNOMAD; peek(GNO)

TARGS.HUM = PGS.DATAS.HUMAN1;
TARGS.MK1 = PGS.DATAS.MOUSE1;
TARGS.MK2 = PGS.DATAS.MOUSE2;



i = sum(GNO.GENE == TARGS.HUM.HL_gene' ,2)>0;
GNO.HUMAN_TARGETS = zeros(height(GNO),1);
GNO.HUMAN_TARGETS(i) = 1;




i = sum(GNO.GENE == TARGS.MK1.ML1_gene' ,2)>0;
GNO.MOUSE1_GENES      = zeros(height(GNO),1);
GNO.MOUSE1_TARGETS    = zeros(height(GNO),1);
GNO.MOUSE1_TARGETS(i) = 1;



i = sum(GNO.GENE == TARGS.MK2.ML2_hs_gene' ,2)>0;
GNO.MOUSE2_GENES      = zeros(height(GNO),1);
GNO.MOUSE2_TARGETS    = zeros(height(GNO),1);
GNO.MOUSE2_TARGETS(i) = 1;





% COUNT NUMBER OF HUMAN TARGETS & MOUSE TARGETS
%------------------------------------------------------------
[C,ia,ic] = unique(GNO.GENE,'stable');
TAB = GNO(ia,:);

clc
fprintf('GNOMAD LoF GENES:	%9.0f \n' ,  height(TAB))
fprintf('~DAWES:            %10.0f \n' , sum(~TAB.DAW_isdawes))
fprintf(' DAWES:            %10.0f \n' , sum(TAB.DAW_isdawes))
fprintf(' HUMAN:            %10.0f \n' , sum(TAB.HUMAN_TARGETS))
fprintf(' MOUSE1:           %10.0f \n' , sum(TAB.MOUSE1_TARGETS))
fprintf(' MOUSE2:           %10.0f \n' , sum(TAB.MOUSE2_TARGETS))
fprintf(' DAWES & HUMAN:    %10.0f \n' , sum(TAB.DAW_isdawes & TAB.HUMAN_TARGETS))
fprintf(' DAWES & MOUSE1:   %10.0f \n' , sum(TAB.DAW_isdawes & TAB.MOUSE1_TARGETS))
fprintf(' DAWES & MOUSE2:   %10.0f \n' , sum(TAB.DAW_isdawes & TAB.MOUSE2_TARGETS))
fprintf(' HUMAN & MOUSE1:   %10.0f \n' , sum(TAB.HUMAN_TARGETS & TAB.MOUSE1_TARGETS))
fprintf(' HUMAN & MOUSE2:   %10.0f \n' , sum(TAB.HUMAN_TARGETS & TAB.MOUSE2_TARGETS))


%------------------------------------------------------------
% sum(TAB.DAW_isdawes    & TAB.DAW_has_mim)
% sum(TAB.HUMAN_TARGETS  & TAB.DAW_has_mim)
% sum(TAB.MOUSE1_TARGETS & TAB.DAW_has_mim)
% sum(TAB.MOUSE2_TARGETS & TAB.DAW_has_mim)
%------------------------------------------------------------







GNOMAD = GNO;






%==========================================================================
%% CLEAN HUMAN DATA TABLE
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD DAW



HUM = PGS.HUM; peek(HUM)
%HUM = convertvars(HUM,@iscell,'string');


% SOME ROWS OF HUM ARE DUPLICATES APART FROM THE HUM.matches ENTRY.
% CONCATINATE THESE ENTRIES FOR ANY DUPLICATE MIM NUMBERS AND REMOVE
% THE DUPLICATE ROW.
%------------------------------------------------------------
[mim,ia,ic] = unique(HUM.mim,'stable');

for i = 1:numel(mim)

    T = HUM(HUM.mim == mim(i) ,:);

    if height(T) > 1

        m = strcat(T.matches(1), "; ", T.matches(2));
    
        T.matches{1} = m;
        T.matches{2} = m;

    end

    HUM(HUM.mim == mim(i) ,:) = T;

end


HUM = HUM(ia,:);
%}





%==========================================================================
%% INSERT DATA FROM HUMAN TABLE INTO GNOMAD TABLE
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD DAW HUM



% ADD COLUMNS TO GNOMAD TABLE TO ACCOMODATE HUM VARIABLES
%------------------------------------------------------------
% HUM(HUM.search_field ~= "tx_clinical_features" ,:) = [];


GNO = GNOMAD;
GNO.mim_number = double(GNO.mim_number);


GNO.HUM_MIM         = zeros(height(GNO),1);
GNO.HUM_GENE        = strings(height(GNO),1);
GNO.HUM_RECESSIVE   = strings(height(GNO),1);
GNO.HUM_AUTOSOMAL   = strings(height(GNO),1);
GNO.HUM_SOMATIC     = strings(height(GNO),1);
GNO.HUM_XLINKED     = strings(height(GNO),1);
GNO.HUM_MOSAIC      = strings(height(GNO),1);
GNO.HUM_DIGENIC     = strings(height(GNO),1);
GNO.HUM_MITOCHON    = strings(height(GNO),1);
GNO.HUM_LETHALITY   = strings(height(GNO),1);






% ADD HUMAN DATA TO GNOMAD TABLE BASED ON MATCHING MIM NUMBER
%------------------------------------------------------------

HUM = convertvars(HUM,@iscell,'string');

[ai,bi] = ismember( GNO.mim_number , HUM.gene_mim_nums );


GNO.HUM_MIM(ai)         = HUM.gene_mim_nums(bi(ai),:);
GNO.HUM_GENE(ai)        = HUM.gene(bi(ai),:);
GNO.HUM_RECESSIVE(ai)   = HUM.recessive(bi(ai),:);
GNO.HUM_AUTOSOMAL(ai)   = HUM.autosomal(bi(ai),:);
GNO.HUM_SOMATIC(ai)     = HUM.somatic(bi(ai),:);
GNO.HUM_XLINKED(ai)     = HUM.xlinked(bi(ai),:);
GNO.HUM_MOSAIC(ai)      = HUM.mosaic(bi(ai),:);
GNO.HUM_DIGENIC(ai)     = HUM.digenic(bi(ai),:);
GNO.HUM_MITOCHON(ai)    = HUM.mito(bi(ai),:);
GNO.HUM_LETHALITY(ai)   = HUM.matches(bi(ai),:);


numel(unique(GNO.HUM_MIM(GNO.HUM_MIM>0)))
%}








%==========================================================================
%% TAG PRENATAL LETHAL GENES
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD DAW HUM GNO



Hmim = HUM.gene_mim_nums(HUM.gene_mim_nums>0);

i = sum(GNO.mim_number == Hmim' ,2)>0; sum(i)

GNO.TF_HAS_HUMAN_DATA = i;

T = regexpi(GNO.HUM_LETHALITY,'embryonic');
L(:,1) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'embryo');
L(:,2) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'prenatal');
L(:,3) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'neonatal');
L(:,4) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'perinatal');
L(:,5) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'fetal');
L(:,6) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'fetus');
L(:,7) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'abortion');
L(:,8) = ~(cellfun(@isempty,T));
T = regexpi(GNO.HUM_LETHALITY,'abort');
L(:,9) = ~(cellfun(@isempty,T));

LETH = sum(L,2)>0;

GNO.TF_PRENATAL_LETHAL = LETH;



sum(GNO.TF_PRENATAL_LETHAL)




GNOMAD = GNO;
%}






%==========================================================================
%% FLAG THE 3435 CANDIDATE DEVELOPMENTAL LETHAL HUMAN GENES
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD DAW



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





% i = sum(DAW.DAWES2_S1.DGENE == DAW.DAWES2_S3.GENES' ,2)>0;
% DAW.DAWES2_S1.CANDIDATE_GENE = i;
%}









%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD



DATAS = PGS.DATAS;

DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP057_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');


