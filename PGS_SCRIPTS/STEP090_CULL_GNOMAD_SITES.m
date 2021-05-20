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
P.fun    = [P.home filesep 'PGS_FUNS'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------










%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP070_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;



GNOMAD.isSV = GNOMAD.VARIANT_CLASS == "SV";



%==========================================================================
%% GET GENE AND VARIANT REPORT
%==========================================================================
clc; clearvars -except P PGS GNOMAD



gene_report(GNOMAD)






%==========================================================================
%% REMOVE SITES WITH FEW SEQUENCED ALLELES (AN)
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)



% max(GNOMAD.AN)
% round(median(GNOMAD.AN))
% round(mean(GNOMAD.AN))
% round(quantile(GNOMAD.AN,.01))
% min(GNOMAD.AN)




close all; histogram(GNOMAD.AN(GNOMAD.allele_type == "SV"));
close all; histogram(GNOMAD.AN(GNOMAD.allele_type ~= "SV"));


GNOMAD( (GNOMAD.AN < 10000)  & (GNOMAD.allele_type == "SV") , :) = [];

GNOMAD( (GNOMAD.AN < 100000) & (GNOMAD.allele_type ~= "SV") , :) = [];








GNOMAD = GENEijk(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% REMOVE SEX CHROMOSOMES
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)




%unique(GNOMAD.CHR)
GNOMAD(GNOMAD.CHR > 22 , :) = [];





GNOMAD = GENEijk(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS THAT HAVE OVER X% MAF
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)




GNOMAD(GNOMAD.AF > .50 , :) = [];





GNOMAD = GENEijk(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS THAT HAVE MORE THAN N% HOMOZYGOUS MINOR ALLELES
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)




pct_nhomalt = (GNOMAD.nhomalt ./ GNOMAD.AC);

i = pct_nhomalt > .10; sum(i)

GNOMAD(i , :) = [];






GNOMAD = GENEijk(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% REMOVE ROWS FLAGGED BY GNOMAD PIPELINE
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)



GNO = GNOMAD;



sum(GNO.isSV)
GNO(GNO.segdup == 1 , :)          = [];
GNO(GNO.lcr == 1 , :)             = [];
GNO(GNO.decoy == 1 , :)           = [];
GNO(GNO.LoF_flags ~= "NaN", :)    = [];
sum(GNO.isSV)



GNOMAD = GNO;


GNOMAD = GENEijk(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% ONLY USE CANONICAL TRANSCRIPTS
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)





sum(GNOMAD.isSV)
GNOMAD( GNOMAD.CANONICAL ~= 1, : ) = [];
sum(GNOMAD.isSV)







GNOMAD = GENEijk(GNOMAD);
gene_report(GNOMAD)
%==========================================================================
%% REMOVE SITES THAT APPEAR LIKE DUPLICATE ALLELES (RUN AS LAST EXCLUSION)
%==========================================================================
clc; clearvars -except P PGS GNOMAD
gene_report(GNOMAD)
sum(GNOMAD.isSV)





[C,ia,ic] = unique([GNOMAD.CHRPOS GNOMAD.AC],'rows','stable');
GNOMAD = GNOMAD(ia,:);






GNOMAD = GENEijk(GNOMAD);
sum(GNOMAD.isSV)
gene_report(GNOMAD)
%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DATAS = PGS.DATAS;
DAWES = PGS.DAWES;


save([P.mat P.f 'PGS_STEP090_OUTPUT.mat'],'GNOMAD','DATAS','DAWES');


gene_report(GNOMAD)









