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
CLINVAR = PGS.CLINVAR;
HUMAN   = PGS.HUMAN;
MOUSE   = PGS.MOUSE;







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
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE


GNOMAD.filter{1253}
GNOMAD.filter{1254}
GNOMAD.filter{1255}
GNOMAD.filter{22637}
GNOMAD.filter{22638}


FIL = GNOMAD.filter;
TER = repmat("NA",height(GNOMAD),1);


% GET CLASS OF EACH CELL
CLASS = cellfun(@class,FIL,'UniformOutput',false);
CLASS = string(CLASS);
gcounts( CLASS )
[Ci,Cnames] = findgroups(CLASS);



% CONVERT ANY CHAR CLASS TO STRING CLASS
TER(Ci == 2) = string(FIL(Ci == 2));
gcounts( TER )



% CONVERT ANY CELL CLASS TO STRING CLASS
[Ci,Cnames] = findgroups(CLASS);
F = FIL(Ci == 1);
% fcell = @(x) strjoin(x);
A = cellfun(@strjoin,F,'UniformOutput',false);
gcounts( A )
B = string(A);
TER(Ci == 1) = B;


gcounts( TER )



GNOMAD.filter = TER;





%==========================================================================
%% GNOMAD.inSegDup | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.inSegDup;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.inSegDup = TF;




%==========================================================================
%% GNOMAD.inLowComplex | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.inLowComplex;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.inLowComplex = TF;





%==========================================================================
%% GNOMAD.inDecoy | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.inDecoy;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.inDecoy = TF;





%==========================================================================
%% GNOMAD.nonpar | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.nonpar;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.nonpar = TF;






%==========================================================================
%% GNOMAD.rf_positive_label | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.rf_positive_label;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.rf_positive_label = TF;






%==========================================================================
%% GNOMAD.rf_negative_label | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.rf_negative_label;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.rf_negative_label = TF;





%==========================================================================
%% GNOMAD.rf_train | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.rf_train;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.rf_train = TF;





%==========================================================================
%% GNOMAD.was_mixed | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.was_mixed;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.was_mixed = TF;




%==========================================================================
%% GNOMAD.has_star | CONVERT STRING "TRUE" OR "FALSE" TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.has_star;

gcounts( VAR )

TF = VAR == "true";

GNOMAD.has_star = TF;




%==========================================================================
%% GNOMAD.strand | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.strand;

gcounts( VAR )

TF = VAR == "1";

GNOMAD.strand = TF;






%==========================================================================
%% GNOMAD.canonical | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.canonical;

gcounts( VAR )

TF = VAR == "YES";

GNOMAD.canonical = TF;






%==========================================================================
%% GNOMAD.genePheno | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE



VAR = GNOMAD.genePheno;

gcounts( VAR )

TF = VAR == "1";

GNOMAD.genePheno = TF;







%==========================================================================
%% GNOMAD.genePheno | CONVERT STRING TO LOGICAL CLASS
%==========================================================================
clc; clearvars -except P GNOMAD CLINVAR HUMAN MOUSE




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
clc; clearvars -except P GNOMAD HUMAN MOUSE CLINVAR

save([P.mat P.f 'PGS_STEP050_OUTPUT.mat'],'GNOMAD','HUMAN','MOUSE','CLINVAR');

writetable(GNOMAD,[P.csv P.f 'PGS_STEP050_OUTPUT.csv'])

