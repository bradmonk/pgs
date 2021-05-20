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
%% GNOMAD BROAD INSTITUTE & MCARTHUR LAB NOTES
%==========================================================================
%{

%------------------------------------------------
FAQ FLAGS
https://gnomad.broadinstitute.org/faq
%------------------------------------------------



Flags that appear on variant pages:
%---------------
AC0: 
    The allele count is zero after filtering out low-confidence genotypes
    (GQ < 20; DP < 10; and AB < 0.2 for het calls) 


AS-VQSR (gnomAD v3 only):
    Failed GATK Allele-Specific Variant Quality Recalibration (AS-VQSR)


InbreedingCoeff: 
    The Inbreeding Coefficient is < -0.3 


RF (gnomAD v2 only):
    Failed random forest filtering thresholds: .055 for exome SNVs, 0.206
    for exome indels, 0.263 genome SNVs, and 0.222 genome indels 




Flags that appear in variant table on gene/region pages:
%---------------
MNV: 
    Multinucleotide variant: the variant is found in phase with another
    variant, altering its functional interpretation 


LCR: 
    Found in a low complexity region: these regions were identified with 
    the symmetric DUST algorithm at a score threshold of 30


LC pLoF:
    LOFTEE Low-confidence pLoF: variant unlikely to be LoF


pLoF Flag: 
    Flagged by LOFTEE; use caution when interpreting the transcript or 
    variant NC Transcript: Marked in a putative LoF category by VEP 
    (essential splice, stop-gained, or frameshift) but appears on a 
    non-protein-coding transcript


%}



%==========================================================================
%% LOAD DATASET
%==========================================================================
clc; clearvars -except P


PGS = load([P.mat P.f 'PGS_STEP040_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;







%==========================================================================
%% ADD VARIABLE DESCRIPTIONS
%==========================================================================
clc; clearvars -except P PGS GNOMAD




GNOMAD.Properties.VariableDescriptions{'AF'} = 'Alternate allele frequency in samples; AF = AC/AN';
GNOMAD.Properties.VariableDescriptions{'VCR'} = 'Variant Carrier Rate; VCR=(AC-NHOMALT)/(.5*AN)';
GNOMAD.Properties.VariableDescriptions{'GCR'} = 'Gene Carrier Rate; GCR = 1 - prod(1 - VCR(i:j));';
GNOMAD.Properties.VariableDescriptions{'GENEi'} = 'Unique gene integer';
GNOMAD.Properties.VariableDescriptions{'GENEj'} = 'Indicates first var in a gene 0/1';
GNOMAD.Properties.VariableDescriptions{'GENEk'} = 'Indicates first var in a gene 1,2,3...';







%==========================================================================
%%                      POS = POS + 1
%==========================================================================
clc; clearvars -except P PGS GNOMAD



GNOMAD.POS = GNOMAD.POS + 1;

GNOMAD.CHRPOS = makeCHRPOS(GNOMAD.CHR,GNOMAD.POS);

GNOMAD = movevars(GNOMAD,{'CHRPOS'},'Before','CHR');








%==========================================================================
%% FIX GNOM.filter, WHICH CONTAINS 1xN CELLS (SHOULD BE STRING)
%==========================================================================
% FOR SOME REASON class(GNOM.filter) = CELL; CHANGE CLASS TO STRING...
% 
% >> GNOM.filter = string(GNOM.filter);
% 
%   Error using string
%   Conversion from cell failed. Element 1254 must be convertible 
%   to a string scalar. 
% 
% 
% >> GNOM.filter{1254}
%   ans =
%       2×1 cell array
%         {'InbreedingCoeff'}
%         {'RF'             }
% 
% 
% THE 'InbreedingCoeff' FILTER FLAG SEEMS ODD. I AM GOING TO REMOVE.
% 
% >> GNOM.filter{1254}(1) = [];
% 
% 
%--------------------------------------------------------------------------
clc; clearvars -except P PGS GNOMAD



unique(GNOMAD.filter)


%gcounts( GNOMAD.filter )
%GNOMAD.filter = TER;





%==========================================================================
%% GNOMAD.inSegDup | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.segdup;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.segdup = VAR;




%==========================================================================
%% GNOMAD.inLowComplex | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.lcr;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.lcr = VAR;





%==========================================================================
%% GNOMAD.inDecoy | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.decoy;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.decoy = VAR;





%==========================================================================
%% GNOMAD.nonpar | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.nonpar;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.nonpar = VAR;






%==========================================================================
%% GNOMAD.rf_positive_label | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.rf_positive_label;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.rf_positive_label = VAR;




%==========================================================================
%% GNOMAD.rf_negative_label | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.rf_negative_label;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.rf_negative_label = VAR;





%==========================================================================
%% GNOMAD.rf_train | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD





VAR = GNOMAD.rf_train;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.rf_train = VAR;






%==========================================================================
%% GNOMAD.was_mixed | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.was_mixed;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.was_mixed = VAR;





%==========================================================================
%% GNOMAD.has_star | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.has_star;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.has_star = VAR;




%==========================================================================
%% GNOMAD.strand | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.STRAND;

gcounts( VAR )

VAR = double(VAR);

gcounts( VAR )

GNOMAD.STRAND = VAR;







%==========================================================================
%% GNOMAD.canonical | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.CANONICAL;

gcounts( VAR )

TF = VAR == "YES";

GNOMAD.CANONICAL = TF;






%==========================================================================
%% GNOMAD.genePheno | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



VAR = GNOMAD.GENE_PHENO;

gcounts( VAR )

TF = VAR == "1";

GNOMAD.GENE_PHENO = TF;




% return
%==========================================================================
%% SORT GNOMAD TABLE AND FILL MISSING GENES WITH PREV
%==========================================================================
clc; clearvars -except P PGS GNOMAD



GNOMAD.CHRPOS = makeCHRPOS(GNOMAD.CHR, GNOMAD.POS);

GNOMAD = movevars(GNOMAD,'CHRPOS','Before','CHR');

GNOMAD = sortrows(GNOMAD,'CHRPOS');

GNOMAD.GENE = fillmissing(GNOMAD.GENE,'previous');

GNOMAD.isLOF = ones(size(GNOMAD,1),1);





%==========================================================================
%% GNOMAD.genePheno | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P PGS GNOMAD




GNOMAD.Properties.VariableDescriptions{'VCR'} = 'Variant Carrier Rate; VCR=(AC-NHOMALT)/(.5*AN)';
GNOMAD.Properties.VariableDescriptions{'GCR'} = 'Gene Carrier Rate; GCR = 1 - prod(1 - VCR(i:j));';
GNOMAD.Properties.VariableDescriptions{'GENEi'} = 'Unique gene integer';


hint(GNOMAD,'VCR');
hint(GNOMAD,'GCR');
hint(GNOMAD,'GENEi');

GNOMAD = hint(GNOMAD, 'VCR'  , 'Variant Carrier Rate; VCR=(AC-NHOMALT)/(.5*AN)');
GNOMAD = hint(GNOMAD, 'GCR'  , 'Gene Carrier Rate; GCR = 1 - prod(1 - VCR(i:j));');
GNOMAD = hint(GNOMAD, 'GENEi', 'Unique gene integer');









%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD

DATAS = PGS.DATAS;


save([P.mat P.f 'PGS_STEP043_OUTPUT.mat'],'GNOMAD','DATAS');


% writetable(GNOMAD,[P.csv P.f 'PGS_STEP050_OUTPUT.csv'])

