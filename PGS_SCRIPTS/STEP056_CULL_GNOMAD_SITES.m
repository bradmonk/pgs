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


PGS = load([P.mat P.f 'PGS_STEP055_B_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;



% GNOMAD DATASET COLUMN REFERENCE
%------------------------------------------------------------
% GNOMAD COLUMNS    T.GENEi            -    T.isSV
% PLI    COLUMNS    T.pLI              -    T.oe_lof
% DAWES  COLUMNS    T.hgnc_id          -    T.isdawes
% MK1    COLUMNS    T.MKO1_TARGETS     -    T.placenta
% MK2    COLUMNS    T.HOMOLO_GENE_ID   -    T.MKO2_TARGETS
% SCR    COLUMNS    T.INVITAE          -    T.MYRIAD
% TARG   COLUMNS    T.HAS_OMIM         -    T.MOUSE2_TARGETS
%------------------------------------------------------------






%==========================================================================
%% REMOVE SITES WITH FEW SEQUENCED ALLELES (AN)
%==========================================================================
clc; clearvars -except P PGS GNOMAD



GNO = GNOMAD;

close all; histogram(GNO.AN)

N = round(quantile(GNO.AN,.005));

GNO( GNO.AN < N , :) = [];



close all; histogram(GNO.AN)


GNOMAD = GNO;





Nreport(GNOMAD)
%==========================================================================
%% REMOVE SEX CHROMOSOMES
%==========================================================================
clc; clearvars -except P PGS GNOMAD



unique(GNOMAD.CHR)


GNOMAD(GNOMAD.CHR == 23 , :) = [];
GNOMAD(GNOMAD.CHR == 25 , :) = [];









Nreport(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS THAT HAVE OVER 50% MAF
%==========================================================================
clc; clearvars -except P PGS GNOMAD




close all; histogram(GNOMAD.AF(GNOMAD.isSV==1 & GNOMAD.AF>.05),.05:.01:1);
hold on;   histogram(GNOMAD.AF(GNOMAD.isLOF==1 & GNOMAD.AF>.05),.05:.01:1);

GNOMAD(GNOMAD.AF > .50 , :) = [];






Nreport(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS THAT HAVE MORE THAN N% HOMOZYGOUS MINOR ALLELES
%==========================================================================
clc; clearvars -except P PGS GNOMAD






NHOMALT_PCT = (GNOMAD.NHOMALT ./ GNOMAD.AC);


close all; histogram( NHOMALT_PCT(GNOMAD.isSV==1  ) , 0:.01:.5);
hold on;   histogram( NHOMALT_PCT(GNOMAD.isLOF==1 ) , 0:.01:.5);

close all; histogram( NHOMALT_PCT(GNOMAD.isSV==1 & NHOMALT_PCT>.01  ) , 0:.01:.5);
hold on;   histogram( NHOMALT_PCT(GNOMAD.isLOF==1 & NHOMALT_PCT>.01 ) , 0:.01:.5);




i = NHOMALT_PCT > .1;

GNOMAD(i , :) = [];







Nreport(GNOMAD)
%==========================================================================
%% REMOVE ROWS FLAGGED BY GNOMAD PIPELINE
%==========================================================================
clc; clearvars -except P PGS GNOMAD


STR = GNOMAD(GNOMAD.isSV==1,:);
LOF = GNOMAD(GNOMAD.isLOF==1,:);


gcounts( STR.inSegDup )
gcounts( STR.inLowComplex )
gcounts( STR.inDecoy )
gcounts( STR.lof_flags ~= "NaN" )

gcounts( LOF.inSegDup )
gcounts( LOF.inLowComplex )
gcounts( LOF.inDecoy )
gcounts( LOF.lof_flags ~= "NaN" )



GNOMAD(GNOMAD.inSegDup == 1 , :)        = [];
GNOMAD(GNOMAD.inLowComplex == 1 , :)    = [];
GNOMAD(GNOMAD.inDecoy == 1 , :)         = [];
GNOMAD(GNOMAD.lof_flags ~= "NaN" & GNOMAD.isLOF==1, :)   = [];





Nreport(GNOMAD)
%==========================================================================
%% REMOVE VARIANTS WITH MAF > N
%==========================================================================
clc; clearvars -except P PGS GNOMAD




close all; histogram(GNOMAD.AF(GNOMAD.isSV==1 & GNOMAD.AF>.02),.02:.01:.5);
hold on;   histogram(GNOMAD.AF(GNOMAD.isLOF==1 & GNOMAD.AF>.02),.02:.01:.5);

GNOMAD(GNOMAD.AF > .30 , :) = [];


close all; histogram(GNOMAD.AF(GNOMAD.isSV==1 & GNOMAD.AF>.02),.02:.01:.5);
hold on;   histogram(GNOMAD.AF(GNOMAD.isLOF==1 & GNOMAD.AF>.02),.02:.01:.5);


Nreport(GNOMAD)
%==========================================================================
%% ONLY USE CANONICAL TRANSCRIPTS
%==========================================================================
% clc; clearvars -except P PGS GNOMAD
% 
% 
% GNOMAD( GNOMAD.canonical ~= 1, : ) = [];
% 
% DMAD( DMAD.canonical ~= 1, : ) = [];
% 
% 
% height(GNOMAD)
% height(DMAD)
% 
% 
%==========================================================================
%% REMOVE SITES THAT APPEAR LIKE DUPLICATE ALLELES (RUN AS LAST EXCLUSION)
%==========================================================================
clc; clearvars -except P PGS GNOMAD




GNOM = GNOMAD;
[C,ia,ic] = unique([GNOM.CHR GNOM.POS GNOM.AC],'rows','stable');
GNOM = GNOM(ia,:);

GNOM = GENEijk(GNOM);



GNOMAD = GNOM;


Nreport(GNOMAD)
%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


GAD = PGS.GNOMAD;

TABLES = PGS.TABLES;

save([P.mat P.f 'PGS_STEP056_OUTPUT.mat'],'GNOMAD','TABLES','GAD');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP056_OUTPUT.csv'])


