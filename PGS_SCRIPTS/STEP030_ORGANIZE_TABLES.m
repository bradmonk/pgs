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
%% LOAD LOFGENES TABLE FROM LOFGENES.mat
%==========================================================================
clearvars -except P


M = load([P.mat P.f 'PGS_STEP020_OUTPUT.mat']);

GNOMAD = M.GNOM;







%==========================================================================
%% LOAD PUBLISHED LETHAL GENE LISTS FROM CSV FILES
%==========================================================================
clearvars -except P GNOMAD


o = detectImportOptions([P.csv '/aleks/PGS_MAIN.xlsx'],'Sheet','HUMAN_LETHAL');
HUMAN1 = readtable([P.csv P.f '/aleks/PGS_MAIN.xlsx'],o);

o = detectImportOptions([P.csv '/aleks/PGS_MAIN.xlsx'],'Sheet','LGLIST_HUMAN');
HUMAN2 = readtable([P.csv P.f '/aleks/PGS_MAIN.xlsx'],o);

o = detectImportOptions([P.csv '/aleks/PGS_MAIN.xlsx'],'Sheet','MOUSE_LETHAL_A');
MOUSE1 = readtable([P.csv P.f '/aleks/PGS_MAIN.xlsx'],o);

o = detectImportOptions([P.csv '/aleks/PGS_MAIN.xlsx'],'Sheet','MOUSE_LETHAL_B');
MOUSE2 = readtable([P.csv P.f '/aleks/PGS_MAIN.xlsx'],o);

o = detectImportOptions([P.csv '/aleks/PGS_MAIN.xlsx'],'Sheet','LGLIST_MOUSE');
MOUSE3 = readtable([P.csv P.f '/aleks/PGS_MAIN.xlsx'],o);

o = detectImportOptions([P.csv '/aleks/PGS_MAIN.xlsx'],'Sheet','CLINVAR');
GUO = readtable([P.csv P.f '/aleks/PGS_MAIN.xlsx'],o);



DATAS.HUMAN1    = HUMAN1;
DATAS.HUMAN2    = HUMAN2;
DATAS.MOUSE1    = MOUSE1;
DATAS.MOUSE2    = MOUSE2;
DATAS.MOUSE3    = MOUSE3;
DATAS.GUO       = GUO;





%==========================================================================
%% CLEAN HUMAN & MOUSE GENE TABLE
%==========================================================================
clc; clearvars -except P GNOMAD DATAS


fprintf('\n\nGNOMAD Variables:\n------------------------\n')
disp(char(string(GNOMAD.Properties.VariableNames')));

fprintf('\n\nHUMAN1 Variables:\n------------------------\n')
disp(char(string(DATAS.HUMAN1.Properties.VariableNames')));

fprintf('\n\nHUMAN2 Variables:\n------------------------\n')
disp(char(string(DATAS.HUMAN2.Properties.VariableNames')));

fprintf('\n\nMOUSE1 Variables:\n------------------------\n')
disp(char(string(DATAS.MOUSE1.Properties.VariableNames')));

fprintf('\n\nMOUSE2 Variables:\n------------------------\n')
disp(char(string(DATAS.MOUSE2.Properties.VariableNames')));

fprintf('\n\nMOUSE3 Variables:\n------------------------\n')
disp(char(string(DATAS.MOUSE3.Properties.VariableNames')));

fprintf('\n\nGUO Variables:\n------------------------\n')
disp(char(string(DATAS.GUO.Properties.VariableNames')));







%==========================================================================
%% RENAME VARIABLES IN GNOMAD TABLE
%==========================================================================
clc; clearvars -except P GNOMAD DATAS


GNOMAD.Properties.VariableNames{'Gene'} = 'GENE_ID';
GNOMAD.Properties.VariableNames{'SYMBOL'} = 'GENE';
GNOMAD.Properties.VariableNames{'start_position'} = 'POS';
GNOMAD.Properties.VariableNames{'end_position'} = 'POS_END';
GNOMAD.Properties.VariableNames{'reference_bases'} = 'REF';
GNOMAD.Properties.VariableNames{'alt'} = 'ALT';
GNOMAD.Properties.VariableNames{'variant_type'} = 'VTYPE';











%==========================================================================
%% REMOVE GENDER COUNTS FROM SUBPOPULATIONS
%==========================================================================
clc; clearvars -except P GNOMAD DATAS



GNOMAD.AN_afr_female           = [];
GNOMAD.AN_afr_male             = [];
GNOMAD.AN_amr_female           = [];
GNOMAD.AN_amr_male             = [];
GNOMAD.AN_asj_female           = [];
GNOMAD.AN_asj_male             = [];
GNOMAD.AN_eas_female           = [];
GNOMAD.AN_eas_male             = [];
GNOMAD.AN_fin_female           = [];
GNOMAD.AN_fin_male             = [];
GNOMAD.AN_nfe_female           = [];
GNOMAD.AN_nfe_male             = [];
GNOMAD.AN_oth_female           = [];
GNOMAD.AN_oth_male             = [];
GNOMAD.AN_sas_female           = [];
GNOMAD.AN_sas_male             = [];
GNOMAD.AC_afr_female           = [];
GNOMAD.AC_afr_male             = [];
GNOMAD.AC_amr_female           = [];
GNOMAD.AC_amr_male             = [];
GNOMAD.AC_asj_female           = [];
GNOMAD.AC_asj_male             = [];
GNOMAD.AC_eas_female           = [];
GNOMAD.AC_eas_male             = [];
GNOMAD.AC_fin_female           = [];
GNOMAD.AC_fin_male             = [];
GNOMAD.AC_nfe_female           = [];
GNOMAD.AC_nfe_male             = [];
GNOMAD.AC_oth_female           = [];
GNOMAD.AC_oth_male             = [];
GNOMAD.AC_sas_female           = [];
GNOMAD.AC_sas_male             = [];
GNOMAD.AF_afr_female           = [];
GNOMAD.AF_afr_male             = [];
GNOMAD.AF_amr_female           = [];
GNOMAD.AF_amr_male             = [];
GNOMAD.AF_asj_female           = [];
GNOMAD.AF_asj_male             = [];
GNOMAD.AF_eas_female           = [];
GNOMAD.AF_eas_male             = [];
GNOMAD.AF_fin_female           = [];
GNOMAD.AF_fin_male             = [];
GNOMAD.AF_nfe_female           = [];
GNOMAD.AF_nfe_male             = [];
GNOMAD.AF_oth_female           = [];
GNOMAD.AF_oth_male             = [];
GNOMAD.AF_sas_female           = [];
GNOMAD.AF_sas_male             = [];

GNOMAD.nhomalt_afr_female      = [];
GNOMAD.nhomalt_afr_male        = [];
GNOMAD.nhomalt_amr_female      = [];
GNOMAD.nhomalt_amr_male        = [];
GNOMAD.nhomalt_asj_female      = [];
GNOMAD.nhomalt_asj_male        = [];
GNOMAD.nhomalt_eas_female      = [];
GNOMAD.nhomalt_eas_male        = [];
GNOMAD.nhomalt_fin_female      = [];
GNOMAD.nhomalt_fin_male        = [];
GNOMAD.nhomalt_nfe_female      = [];
GNOMAD.nhomalt_nfe_male        = [];
GNOMAD.nhomalt_oth_female      = [];
GNOMAD.nhomalt_oth_male        = [];
GNOMAD.nhomalt_sas_female      = [];
GNOMAD.nhomalt_sas_male        = [];



%==========================================================================
%% REORDER VARIABLES IN GNOMAD TABLE
%==========================================================================
clc; clearvars -except P GNOMAD DATAS

peek(GNOMAD)


GNOMAD = movevars(GNOMAD,{'GENE'},'After','ALT');


GNOMAD = movevars(GNOMAD,{'AN_raw','AN_popmax','AN_female','AN_male'},'Before','AN_afr');
GNOMAD = movevars(GNOMAD,{'AC_raw','AC_popmax','AC_female','AC_male'},'Before','AC_afr');
GNOMAD = movevars(GNOMAD,{'AF_raw','AF_popmax','AF_female','AF_male'},'Before','AF_afr');
GNOMAD = movevars(GNOMAD,{'nhomalt_raw','nhomalt_popmax','nhomalt_female','nhomalt_male'},'Before','AF_afr');

GNOMAD = movevars(GNOMAD,{'AA_MAF','AFR_MAF','AMR_MAF','EA_MAF','EAS_MAF','EUR_MAF','SAS_MAF','GMAF'},'After','nhomalt_sas');

peek(GNOMAD)



%==========================================================================
%% CHANGE VARIABLE CLASS
%==========================================================================
clc; clearvars -except P GNOMAD DATAS


% COLUMNS 8-117 SHOULD BE NUMERIC



V = GNOMAD.Properties.VariableNames;

for i = 8:117
    
    if ~isnumeric(GNOMAD.(V{i}))
        GNOMAD.(V{i}) = double(GNOMAD{:,i});
    end
    
end






%==========================================================================
%% RENAME VARIABLES IN GNOMAD TABLE
%==========================================================================
%{
clc; clearvars -except P GNOMAD HUMAN MOUSE CLINVAR
disp(char(string(GNOMAD.Properties.VariableNames')));

T = GNOMAD;


T.Properties.VariableNames{'names'}                 = 'rs';
T.Properties.VariableNames{'BASES'}                 = 'refBase';
T.Properties.VariableNames{'QD'}                    = 'callConf';
T.Properties.VariableNames{'SOR'}                   = 'strandBias';
T.Properties.VariableNames{'MQRankSum'}             = 'Z_mapQuality';
T.Properties.VariableNames{'ReadPosRankSum'}        = 'Z_posBias';
T.Properties.VariableNames{'BaseQRankSum'}          = 'Z_baseQuality';
T.Properties.VariableNames{'ClippingRankSum'}       = 'Z_hardClip';
T.Properties.VariableNames{'VQSLOD'}                = 'oddsTrueVar';
T.Properties.VariableNames{'segdup'}                = 'inSegDup';
T.Properties.VariableNames{'lcr'}                   = 'inLowComplex';
T.Properties.VariableNames{'decoy'}                 = 'inDecoy';
T.Properties.VariableNames{'ALT'}                   = 'altBase';
T.Properties.VariableNames{'variant_type'}          = 'varType';
T.Properties.VariableNames{'GENE'}                  = 'geneID';
T.Properties.VariableNames{'SYMBOL'}                = 'geneName';
T.Properties.VariableNames{'N_ALT_ALLELES'}         = 'N_alleleTypes';
T.Properties.VariableNames{'ALLELE'}                = 'altAllele';
T.Properties.VariableNames{'ALLELE_TYPE'}           = 'alleleType';
T.Properties.VariableNames{'CONSEQUENCE'}           = 'consequence';
T.Properties.VariableNames{'IMPACT'}                = 'impact';
T.Properties.VariableNames{'BIOTYPE'}               = 'biotype';
T.Properties.VariableNames{'InbreedingCoeff'}       = 'inbredCoef';
T.Properties.VariableNames{'DP'}                    = 'depth';
T.Properties.VariableNames{'PAB_MAX'}               = 'pabMax';
T.Properties.VariableNames{'POPMAX'}                = 'popMax';
T.Properties.VariableNames{'AC_POPMAX'}             = 'popMax_AC';
T.Properties.VariableNames{'AN_POPMAX'}             = 'popMax_AN';
T.Properties.VariableNames{'AF_POPMAX'}             = 'popMax_AF';
T.Properties.VariableNames{'NHOMALT_POPMAX'}        = 'popMax_NHOMALT';
T.Properties.VariableNames{'FEATURE_TYPE'}          = 'featureType';
T.Properties.VariableNames{'FEATURE'}               = 'feature';
T.Properties.VariableNames{'EXON'}                  = 'exon';
T.Properties.VariableNames{'INTRON'}                = 'intron';
T.Properties.VariableNames{'PROTEIN_POS'}           = 'proteinPos';
T.Properties.VariableNames{'AMINO_ACIDS'}           = 'aminoAcids';
T.Properties.VariableNames{'CODONS'}                = 'codons';
T.Properties.VariableNames{'EXISTING_VAR'}          = 'existingVar';
T.Properties.VariableNames{'ALLELE_NUM'}            = 'alleleNum';
T.Properties.VariableNames{'STRAND'}                = 'strand';
T.Properties.VariableNames{'FLAGS'}                 = 'flag';
T.Properties.VariableNames{'VARIANT_CLASS'}         = 'varClass';
T.Properties.VariableNames{'MINIMISED'}             = 'minimised';
T.Properties.VariableNames{'SYMBOL_SOURCE'}         = 'symbolSource';
T.Properties.VariableNames{'CANONICAL'}             = 'canonical';
T.Properties.VariableNames{'GENE_PHENO'}            = 'genePheno';
T.Properties.VariableNames{'HGVS_OFFSET'}           = 'HGVS_offset';
T.Properties.VariableNames{'SOMATIC'}               = 'somatic';
T.Properties.VariableNames{'PHENO'}                 = 'pheno';
T.Properties.VariableNames{'LOF'}                   = 'lof';
T.Properties.VariableNames{'LOF_FLAGS'}             = 'lof_flags';



GNOMAD = T;
%}




%==========================================================================
%% REORDER VARIABLES IN GNOMAD TABLE
%==========================================================================
%{
clc; clearvars -except P GNOMAD HUMAN MOUSE CLINVAR
disp(char(string(GNOMAD.Properties.VariableNames')));

T = GNOMAD;


ORDER = [
    1
    2
    3
    60
    61
    9
    4
    36
    57
    5
    6
    7
    8
    37
    38
    40
    42
    43
    44
    45
    46
    47
    48
    49
    50
    51
    10
    39
    58
    59
    64
    11
    12
    13
    14
    15
    16
    17
    18
    19
    20
    21
    22
    23
    24
    25
    26
    27
    28
    29
    30
    31
    32
    33
    34
    35
    41
    52
    53
    54
    55
    56
    62
    63
    65
    66
    67
    68
    69
    70
    71
    72
    73
    74
    75
    76
    77
    78
    79
    80
    81
    82
    83
    84
    85
    86
    87
    88
    89
    90
];


T = T(:,ORDER);



GNOMAD = T;
%}









%==========================================================================
%% ADD VARIABLE DESCRIPTIONS
%==========================================================================
%{



VariableDescriptions = {...
'Chromosome Reference name.';
'Start position (0-based). Corresponds to the first base of the string of reference bases.';
'End position (0-based). Corresponds to the first base after the last base in the reference allele.';
'The gene symbol';
'Ensembl stable ID of affected gene';
'Variant names (e.g. RefSNP ID).';
'Reference bases.';
'Alternate base.';
'The ALT part of the annotation field.';
'Total number of alleles in samples';
'Total number of alleles in samples, before removing low-confidence genotypes';
'Total number of alleles in male samples';
'Total number of alleles in female samples';
'Alternate allele count for samples';
'Alternate allele frequency in samples; AF = AC/AN';
'Total number of alternate alleles observed at variant locus';
'Alternate allele count for samples, before removing low-confidence genotypes';
'Alternate allele frequency in samples, before removing low-confidence genotypes';
'Count of homozygous individuals in samples, before removing low-confidence genotypes';
'Count of homozygous individuals in samples';
'Alternate allele count for male samples';
'Alternate allele count for female samples';
'Alternate allele frequency in male samples';
'Alternate allele frequency in female samples';
'Count of homozygous individuals in male samples';
'Count of homozygous individuals in female samples';
'Variant type (snv, indel, multi-snv, multi-indel, or mixed)';
'Allele type (snv, ins, del, or mixed)';
'Consequence type of this variant';
'The impact modifier for the consequence type';
'Biotype of transcript or regulatory feature';
'Phred-scaled quality score (-10log10 prob(call is wrong)). Higher values imply better quality.';
'List of failed filters (if any) or "PASS" indicating the variant has passed all filters.';
'Random forest prediction probability for a site being a true variant';
'Phred-scaled p-value of Fishers exact test for strand bias';
'Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation';
'Root mean square of the mapping quality of reads across all samples';
'Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities';
'Variant call confidence normalized by depth of sample reads supporting a variant';
'Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias';
'Strand bias estimated by the symmetric odds ratio test';
'Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities';
'Z-score from Wilcoxon rank sum test of alternate vs. reference number of hard clipped bases';
'Depth of informative coverage for each sample; reads with MQ=255 or with bad mates are filtered';
'Log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model';
'Worst-performing annotation in the VQSR Gaussian mixture model';
'Variant falls within a segmental duplication region';
'Variant falls within a low complexity region';
'Variant falls within a reference decoy region';
'Variant (on sex chromosome) falls outside a pseudoautosomal region';
'Variant was labelled as a positive example for training of random forest model';
'Variant was labelled as a negative example for training of random forest model';
'Random forest training label';
'Variant was used in training random forest model';
'Variant type was mixed';
'Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset';
'Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of AB=0.5';
'Population with maximum AF';
'Allele count in the population with the maximum AF';
'Total number of alleles in the population with the maximum AF';
'Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry)';
'Count of homozygous individuals in the population with the maximum allele frequency';
'Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature.';
'Ensembl stable ID of feature';
'The exon number (out of total number)';
'The intron number (out of total number)';
'Relative position of amino acid in protein';
'Reference and variant amino acids. Only given if the variant affects the protein-coding sequence';
'The alternative codons with the variant base in upper case';
'Known identifier of existing variant';
'Allele number from input; 0 is reference, 1 is first alternate etc';
'The DNA strand (1 or -1) on which the transcript/feature lies';
'Transcript quality flags (cds_start_NF, cds_start_NF)';
'Sequence Ontology variant class';
'Alleles in this variant have been converted to minimal representation before consequence calculation';
'The source of the gene symbol';
'HUGO Gene Nomenclature Committee approved symbol';
'A flag indicating if the transcript is denoted as the canonical transcript for this gene';
'The CCDS identifer for this transcript, where applicable';
'The Ensembl protein identifier of the affected transcript';
'Best match UniProtKB/Swiss-Prot accession of protein product';
'Best match UniProtKB/TrEMBL accession of protein product';
'Best match UniParc accession of protein product';
'Indicates if overlapped gene is associated with a phenotype, disease or trait';
'Indicates by how many bases the HGVS notations for this variant have been shifted';
'ClinVar clinical significance of the dbSNP variant';
'Somatic status of existing variant(s); multiple values correspond to multiple values in the Existing_variation field';
'Indicates if existing variant is associated with a phenotype, disease or trait; multiple values correspond to multiple values in the Existing_variation field';
'HC: HIGH CONFIDENCE; LC: LOW CONFIDENCE';
'Loss of Function FLAGS'
};




GNOMAD.Properties.VariableDescriptions = VariableDescriptions;


GNOMAD.Properties.VariableDescriptions
%}




%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P GNOMAD DATAS



save([P.mat P.f 'PGS_STEP030_OUTPUT.mat'],'GNOMAD','DATAS');

% writetable(GNOMAD,[P.csv P.f 'PGS_STEP3_OUTPUT.csv'])







