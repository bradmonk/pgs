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


PGS = load([P.mat P.f 'PGS_STEP057_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;







%==========================================================================
%% CLEAN HUMAN DATA TABLE
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD



HUM = PGS.TARGS.HUM; peek(HUM)
%HUM = convertvars(HUM,@iscell,'string');


% SOME ROWS OF HUM ARE DUPLICATES APART FROM THE HUM.matches ENTRY.
% CONCATINATE THESE ENTRIES FOR ANY DUPLICATE MIM NUMBERS AND REMOVE
% THE DUPLICATE ROW.
%------------------------------------------------------------
% [mim,ia,ic] = unique(HUM.mim,'stable');

[mim,ia,ic] = unique(HUM.HL_mim_num,'stable');

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
%% DAWES DATASET COLUMN DEFINITIONS
%==========================================================================
%{
% DAWES HUMAN CRITERIA
%----------------------------------
% hgnc_id               Gene IDs	HUGO Gene Nomenclature Committee
% gene                  Gene symbols	HUGO Gene Nomenclature Committee
% gene_name             Full Gene names	HUGO Gene Nomenclature Committee
% gene_family           Gene families	HUGO Gene Nomenclature Committee
% *omim                 Y is a clinically relevant OMIM gene (OMIM)
% mim_number            MIM number if the gene is a clinically relevant OMIM gene (OMIM)
% phenotype             Phenotypes associated with the gene in OMIM (OMIM)
% Inheritance_pattern	Inheritance patterns associated with the gene in OMIM (OMIM)
% *human_lethal_B       Y if gene appears in our curated list of human lethal genes (OMIM)
% human_lethal_A        Y if gene appears in our stringently curated list of human lethal genes (OMIM)
% lethal_phen           Lethal phenotypes associated with the gene (OMIM)
% *lethal_inheritance	Inheritance patterns of lethal phenotypes associated with the gene (OMIM)
% gnomAD                Y if the gene has constraint scores (gnomad 2.1)
% mis_z                 missense constraint scores (gnomad 2.1)
% syn_z                 synonymous constraint scores (gnomad 2.1)
% pLI                   pLI scores (gnomad 2.1)
% constrained           Y if the gene has either missense or LoF constraint (gnomad 2.1)
% 
% 
% 
% DAWES MOUSE CRITERIA
%----------------------------------
% *omim                     Y is a clinically relevant OMIM gene (OMIM)
% MGI_ID                    ID of gene in MGI
% lethal_MGI                Y if gene has any lethal MP terms in hom KOs (MGI)
% lethal_MP_ID              IDs of the lethal MP terms for the gene (MGI)
% lethal_MP_phen            Desc of the lethal MP terms for the gene (MGI)
% all_MP_ID                 IDs of all MP terms noted to the gene (MGI)
% all_MP_phen               Desc of all the MP terms noted to the gene (MGI)
% high_MP_ID                IDs of high-level descriptive MP terms for gene (MGI)
% high_MP_phen              Desc of high-level descriptive MP terms for gene (MGI)
% allele_info               Info on type of allele present in hom mouse KO (MGI)
% mouse_symbol              Gene symbol orthologous to human gene symbol (MGI)
% lethal_IMPC               Y if gene has any lethal MP terms in hom KOs (IMPC)
% IMPC_all_MP_ID            IDs of all MP terms for gene (IMPC)
% IMPC_all_MP_phen          Desc of all the MP terms for gene (IMPC)
% IMPC_lethal_MP_ID         IDs of the lethal MP terms for gene (IMPC)
% IMPC_lethal_MP_phen       Desc of the lethal MP terms for gene (IMPC)
% IMPC_ko                   Y if gene has a hom mouse KO (IMPC)
% *mouse_ko                  Y if gene has a hom mouse KO (MGI/IMPC)
% *lethal_mouse              Y if gene is lethal in hom KO mice  (MGI/IMPC)
% lethal_het_MGI            Y if gene has any lethal MP terms in het KOs (MGI)
% lethal_het_MP_ID          IDs of lethal MP terms for gene in MGI het KOs (MGI)
% lethal_het_MP_phen        Desc of lethal MP terms for gene in MGI het KOs (MGI)
% all_het_MP_ID             IDs of all MP terms for gene in MGI het KOs (MGI)
% all_het_MP_phen           Desc of all the MP terms for gene in MGI het Kos (MGI)
% het_allele_info           Info on type of allele present in het KOs (MGI)
% lethal_het_IMPC           Y if gene has any lethal MP terms in het KOs (IMPC)
% IMPC_het_all_MP_ID        IDs of all MP terms for gene in IMPC het KOs (IMPC)
% IMPC_het_all_MP_phen      Desc of all MP terms for gene in IMPC het Kos (IMPC)
% IMPC_het_lethal_MP_phen	Desc of lethal MP terms for gene in IMPC het KOs (IMPC)
% IMPC_het_ko               Y if gene has het KO phenotype data (IMPC)
% mouse_het_ko              Y if gene has het KO phenotype data (MGI/IMPC)
% lethal_het_mouse          Y if gene has any lethal MP terms in het KOs (MGI/IMPC)
%}





%==========================================================================
%% DETERMINE HUMAN & MOUSE TARGET GENES
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DAW = PGS.DAWES.DAWES2_S1;



% HUMAN GENE TARGETS
%----------------------------------
HLT = DAW(DAW.omim == "Y",:);                                % 3187 GENES
HLT = HLT(HLT.human_lethal_B == "Y" ,:);                     %  624 GENES
HLT = HLT(contains(HLT.lethal_inheritance,["AR","XLr"]) ,:); %  504 GENES


% MOUSE GENE TARGETS
%----------------------------------
MLT = DAW(DAW.mouse_ko == "Y",:);                            % 9397 GENES
MLT = MLT(MLT.lethal_mouse == "Y" ,:);                       % 3684 GENES
MLT = MLT(MLT.omim ~= "Y" ,:);                               % 2377 GENES






%==========================================================================
%% DETERMINE HUMAN & MOUSE TARGET GENES WITH GNOMAD LOF MUTATIONS
%==========================================================================
clc; clearvars -except P PGS GNOMAD DAW HLT MLT




% GNOMAD HUMAN GENE TARGETS
%----------------------------------
GHT = GNOMAD(GNOMAD.DAW_has_mim == 1,:);    
sum(GHT.GENEj)                                      % 1280 GENES, 7474 SNPs

GHT = GHT(GHT.DAW_human_lethal_B == 1 ,:);
sum(GHT.GENEj)                                      %  252 GENES, 1118 SNPs

GHT = GHT(contains(GHT.DAW_lethal_inheritance,["AR","XLr"]) ,:);
sum(GHT.GENEj)                                      %  224 GENES, 1001 SNPs






% GNOMAD MOUSE GENE TARGETS
%----------------------------------
GMT = GNOMAD(GNOMAD.DAW_mouse_ko == 1,:);    
sum(GMT.GENEj)                                      % 3581 GENES, 20529 SNPs

GMT = GMT(GMT.DAW_lethal_mouse == 1 ,:);
sum(GMT.GENEj)                                      % 1182 GENES,  6151 SNPs

GMT = GMT(GMT.DAW_has_mim == 0,:);    
sum(GMT.GENEj)                                      %  702 GENES,  3545 SNPs




%==========================================================================
%% MARK TARGETS
%==========================================================================
clc; clearvars -except P PGS GNOMAD DAW HLT MLT GHT GMT


GNO = GNOMAD;

GNO.HUMAN_TARGS = zeros(height(GNO),1);
GNO.MOUSE_TARGS = zeros(height(GNO),1);




% MARK HUMAN TARGETS IN GNOMAD TABLE
%---------------------------------------
i = ismember( GNO.GENE , GHT.GENE );
GNO.HUMAN_TARGS(i) = 1;

sum(GNO.HUMAN_TARGS)                      % 1001 SNPs
sum(GNO.GENEj(GNO.HUMAN_TARGS==1))        %  224 GENES




% MARK MOUSE TARGETS IN GNOMAD TABLE
%---------------------------------------
i = ismember( GNO.GENE , GMT.GENE );
GNO.MOUSE_TARGS(i) = 1;

sum(GNO.MOUSE_TARGS)                      % 3545 SNPs
sum(GNO.GENEj(GNO.MOUSE_TARGS==1))        %  702 GENES




G = GNO(GNO.MOUSE_TARGS==1 | GNO.MOUSE_TARGS==1 ,:);
height(G)                                   % 4546 SNPs  (1001 + 3545)
sum(G.GENEj)                                %  926 GENES (702 + 224)




GNOMAD = GNO;










%==========================================================================
%% INSERT DATA FROM HUMAN TABLE INTO GNOMAD TABLE
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD



% ADD COLUMNS TO GNOMAD TABLE TO ACCOMODATE HUM VARIABLES
%------------------------------------------------------------
% HUM(HUM.search_field ~= "tx_clinical_features" ,:) = [];


GNO = GNOMAD;

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
% HUM(HUM.search_field ~= "tx_clinical_features" ,:) = [];


HUM = PGS.DATAS.HUMAN2;

HUM = convertvars(HUM,@iscell,'string');

[ai,bi] = ismember( GNO.DAW_mim , HUM.gene_mim_nums );


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




GNOMAD = GNO;
%}






%==========================================================================
%% TAG RECESSIVE GENES
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD



GNO = GNOMAD;

HUM = PGS.HUM_OLD;



Hmim = HUM.gene_mim_nums(HUM.gene_mim_nums>0);

i = sum(GNO.DAW_mim == Hmim' ,2)>0; sum(i)

GNO.HUM_has_hum_data = i;

T = regexpi(GNO.HUM_LETHALITY,'recessive');
L(:,1) = ~(cellfun(@isempty,T));

RECESSIVE = sum(L,2)>0;

sum(RECESSIVE)

GNO.HUM_IS_RECESSIVE = RECESSIVE;



sum(GNO.HUM_IS_RECESSIVE)




GNOMAD = GNO;










%==========================================================================
%% TAG PRENATAL LETHAL GENES
%==========================================================================
%{.
clc; clearvars -except P PGS GNOMAD



GNO = GNOMAD;

HUM = PGS.HUM_OLD;



Hmim = HUM.gene_mim_nums(HUM.gene_mim_nums>0);

i = sum(GNO.DAW_mim == Hmim' ,2)>0; sum(i)

GNO.HUM_has_hum_data = i;

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
T = regexpi(GNO.HUM_LETHALITY,'stillborn');
L(:,10) = ~(cellfun(@isempty,T));


LETH = sum(L,2)>0;

GNO.HUM_IS_PRE_LETHAL = LETH;



sum(GNO.HUM_IS_PRE_LETHAL)




GNOMAD = GNO;
%}







%==========================================================================
%% FLAG THE 3435 CANDIDATE DEVELOPMENTAL LETHAL HUMAN GENES
%==========================================================================
%{
clc; clearvars -except P PGS GNOMAD



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





i = sum(PGS.DAWES.DAWES2_S1.DGENE == PGS.DAWES.DAWES2_S3.DGENES' ,2)>0;
PGS.DAWES.DAWES2_S1.CANDIDATE_GENE = i;






%}







%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DATAS = PGS.DATAS;

DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP058_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');





