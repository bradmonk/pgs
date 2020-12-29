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
clearvars -except P


M = load([P.mat P.f 'PGS_STEP030_OUTPUT.mat']);

GNOMAD  = M.GNOMAD;
HUMAN   = M.HUMAN;
MOUSE   = M.MOUSE;
CLINVAR = M.CLINVAR;


unique(GNOMAD.CHR)



%==========================================================================
%% COMPUTE VARIANT CARRIER RATE (VCR)
%==========================================================================
% 
% VCR = (AC - HOM) ./ (.5 .* AN);
% 
% 
% Where...
% AC is the total allele count for the variant, HOM is the number of
% individuals who are homozygous for the variant, and AN is the total
% number of alleles analyzed for the variant.
% 
% AN:  count of all alleles measured (REF + ALT)
% AC:  count of ALT alleles (ALT_HOM + ALT_HET)
% HOM: count of homozygous ALT people (ALT_HOM)
%--------------------------------------------------------------------------
clc; clearvars -except P GNOMAD HUMAN MOUSE CLINVAR




% COMPUTE VCR FOR EACH VARIANT
GNOMAD.VCR = (GNOMAD.AC - GNOMAD.NHOMALT) ./ (.5 .* GNOMAD.AN);







%==========================================================================
%% COMPUTE GENE CARRIER RATE (GCR)
%==========================================================================
% 
% GCR(g) = 1 - prod(1 - VCR(1:k));
% 
% 
% The equation above computes the gene carrier rate (GCR) for a gene 'g',
% where VCR is the variant carrier rates. The product is computed for all 
% variants 1:k of a given gene.
% 
%--------------------------------------------------------------------------
clc; clearvars -except P GNOMAD HUMAN MOUSE CLINVAR




% ADD GCR COLUMN TO TABLE
GNOMAD.GCR = nan(height(GNOMAD),1);




% GET ARRAY INDEX LOCATIONS FOR EACH GENE
[~,~,GNOMAD.GENEi] = unique(GNOMAD.geneID,'stable');





% COMPUTE GCR FOR EACH GENE
for g = 1:numel(unique(GNOMAD.GENEi))

    GNOMAD.GCR(GNOMAD.GENEi==g) = 1 - prod(1 - GNOMAD.VCR(GNOMAD.GENEi==g));

end







%==========================================================================
%% COMPUTE CUMULATIVE CARRIER RATE (CCR)
%==========================================================================
% 
% CCR = 1 - prod(1 - GCR(L) )
% 
% 
% The equation above cumulative carrier (CCR) for a list of genes 'L',
% where GCR is an array of gene carrier rates. The product is computed 
% for the 1-GCR values of all genes in list L.
% 
%--------------------------------------------------------------------------
% % EXAMPLE...
% 
% 
% n = height(GEN);
% 
% GCR = GEN.GCR;
% 
% 
% cc = 1000;
% gg = 30;
% 
% CCR = nan(gg,cc);
% 
% for c = 1:cc
% for g = 1:gg
% 
%     r = randperm(n,g);
% 
%     CCR(g,c) = 1 - prod(1 - GCR(r));
% 
% end
% end
% 
% 
%--------------------------------------------------------------------------







%==========================================================================
%% ADD VARIABLE DESCRIPTIONS
%==========================================================================


GNOMAD.Properties.VariableDescriptions{'AF'} = 'Alternate allele frequency in samples; AF = AC/AN';
GNOMAD.Properties.VariableDescriptions{'VCR'} = 'Variant Carrier Rate; VCR=(AC-NHOMALT)/(.5*AN)';
GNOMAD.Properties.VariableDescriptions{'GCR'} = 'Gene Carrier Rate; GCR = 1 - prod(1 - VCR(i:j));';
GNOMAD.Properties.VariableDescriptions{'GENEi'} = 'Unique gene integer';






%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P GNOMAD HUMAN MOUSE CLINVAR

save([P.mat P.f 'PGS_STEP040_OUTPUT.mat'],'GNOMAD','HUMAN','MOUSE','CLINVAR');

writetable(GNOMAD,[P.csv P.f 'PGS_STEP040_OUTPUT.csv'])














