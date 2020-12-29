%==========================================================================
%% SCRIPT_NAME
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home  = '/Users/bradleymonk/Documents/MATLAB/UCSF/PGS/GNOMAD2';
cd(P.home);
P.srip   = [P.home filesep 'PGS_SCRIPTS'];
P.json   = [P.home filesep 'PGS_JSON'];
P.gnomat = [P.home filesep 'PGS_GNOMAT'];
P.funs   = [P.home filesep 'PGS_FUNS'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;
%------------------------------------------------




%==========================================================================
%% GET FILE PATHS
%==========================================================================
clearvars -except P


P.CHR = 'CHR01';


JSON = jsondecode(fileread([P.json P.f P.CHR '.json']));





%==========================================================================
%% GET FILE PATHS
%==========================================================================
clearvars -except P JSON



TAB = struct2table(JSON);
peek(TAB)


if strcmp(TAB.CHR{1},'X')
    TAB.CHR(:) = {'23'};
end
if strcmp(TAB.CHR{1},'Y')
    TAB.CHR(:) = {'24'};
end



TAB.CHR 					= double(string(char(TAB.CHR)));
TAB.POS 					= double(string(char(TAB.POS)));
TAB.POS_END 				= double(string(char(TAB.POS_END)));
TAB.BASES 					= string(TAB.BASES);
TAB.AN 					    = double(string(char(TAB.AN)));
TAB.AN_raw 					= double(string(char(TAB.AN_raw)));
TAB.AN_male 				= double(string(char(TAB.AN_male)));
TAB.AN_female 				= double(string(char(TAB.AN_female)));
TAB.variant_type 			= string(TAB.variant_type);
TAB.quality 				= double(string(char(TAB.quality)));
TAB.rf_tp                   = double(string(char(TAB.rf_tp)));
TAB.FS 					    = double(string(char(TAB.FS)));
TAB.InbreedingCoeff 		= double(string(char(TAB.InbreedingCoeff)));
TAB.MQ 					    = double(string(char(TAB.MQ)));
TAB.SOR 					= double(string(char(TAB.SOR)));
TAB.DP 					    = double(string(char(TAB.DP)));
TAB.ALT 					= string(TAB.ALT);
TAB.AC 					    = double(string(char(TAB.AC)));
TAB.AF 					    = double(string(char(TAB.AF)));
TAB.ALLELE_TYPE 			= string(TAB.ALLELE_TYPE);
TAB.N_ALT_ALLELES 			= double(string(char(TAB.N_ALT_ALLELES)));
TAB.AC_RAW 					= double(string(char(TAB.AC_RAW)));
TAB.AF_RAW 					= double(string(char(TAB.AF_RAW)));
TAB.NHOMALT_RAW 			= double(string(char(TAB.NHOMALT_RAW)));
TAB.NHOMALT 				= double(string(char(TAB.NHOMALT)));
TAB.AC_MALE 				= double(string(char(TAB.AC_MALE)));
TAB.AC_FEMALE 				= double(string(char(TAB.AC_FEMALE)));
TAB.ALLELE 					= string(TAB.ALLELE);
TAB.CONSEQUENCE 			= string(TAB.CONSEQUENCE);
TAB.IMPACT 					= string(TAB.IMPACT);
TAB.SYMBOL 					= string(TAB.SYMBOL);
TAB.GENE 					= string(TAB.GENE);
TAB.FEATURE_TYPE 			= string(TAB.FEATURE_TYPE);
TAB.FEATURE 				= string(TAB.FEATURE);
TAB.BIOTYPE 				= string(TAB.BIOTYPE);
TAB.ALLELE_NUM 				= double(string(char(TAB.ALLELE_NUM)));
TAB.STRAND 					= string(TAB.STRAND);
TAB.VARIANT_CLASS 			= string(TAB.VARIANT_CLASS);
TAB.MINIMISED 				= double(string(char(TAB.MINIMISED)));
TAB.SYMBOL_SOURCE 			= string(TAB.SYMBOL_SOURCE);
TAB.ENSP 					= string(TAB.ENSP);
TAB.UNIPARC 				= string(TAB.UNIPARC);
TAB.LOF 					= string(TAB.LOF);



%%


EM = cellfun(@isempty,TAB.AF_MALE);              TAB.AF_MALE(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.AF_FEMALE);            TAB.AF_FEMALE(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.NHOMALT_MALE);         TAB.NHOMALT_MALE(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.NHOMALT_FEMALE);       TAB.NHOMALT_FEMALE(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.ClippingRankSum);      TAB.ClippingRankSum(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.BaseQRankSum);         TAB.BaseQRankSum(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.ReadPosRankSum);       TAB.ReadPosRankSum(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.QD);                   TAB.QD(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.MQRankSum);            TAB.MQRankSum(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.PAB_MAX);              TAB.PAB_MAX(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.names);                TAB.names(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.filter);               TAB.filter(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.VQSLOD);               TAB.VQSLOD(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.VQSR_culprit);         TAB.VQSR_culprit(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.segdup);               TAB.segdup(EM) = {NaN};
EM = cellfun(@isempty,TAB.lcr);                  TAB.lcr(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.decoy);                TAB.decoy(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.nonpar);               TAB.nonpar(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.rf_positive_label);    TAB.rf_positive_label(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.rf_negative_label);    TAB.rf_negative_label(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.rf_label);             TAB.rf_label(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.rf_train);             TAB.rf_train(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.was_mixed);            TAB.was_mixed(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.has_star);             TAB.has_star(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.POPMAX);               TAB.POPMAX(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.AC_POPMAX);            TAB.AC_POPMAX(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.AN_POPMAX);            TAB.AN_POPMAX(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.AF_POPMAX);            TAB.AF_POPMAX(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.NHOMALT_POPMAX);       TAB.NHOMALT_POPMAX(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.EXON);                 TAB.EXON(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.INTRON);               TAB.INTRON(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.PROTEIN_POS);          TAB.PROTEIN_POS(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.AMINO_ACIDS);          TAB.AMINO_ACIDS(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.CODONS);               TAB.CODONS(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.EXISTING_VAR);         TAB.EXISTING_VAR(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.FLAGS);                TAB.FLAGS(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.HGNC_ID);              TAB.HGNC_ID(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.CANONICAL);            TAB.CANONICAL(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.CCDS);                 TAB.CCDS(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.SWISSPROT);            TAB.SWISSPROT(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.TREMBL);               TAB.TREMBL(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.GENE_PHENO);           TAB.GENE_PHENO(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.HGVS_OFFSET);          TAB.HGVS_OFFSET(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.CLIN_SIG);             TAB.CLIN_SIG(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.SOMATIC);              TAB.SOMATIC(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.PHENO);                TAB.PHENO(EM) = {'NaN'};
EM = cellfun(@isempty,TAB.LOF_FLAGS);            TAB.LOF_FLAGS(EM) = {'NaN'};





%%

TAB.AF_MALE 				= double(string(char(TAB.AF_MALE)));
TAB.AF_FEMALE 				= double(string(char(TAB.AF_FEMALE)));
TAB.NHOMALT_MALE 			= double(string(char(TAB.NHOMALT_MALE)));
TAB.NHOMALT_FEMALE 			= double(string(char(TAB.NHOMALT_FEMALE)));
TAB.ClippingRankSum 		= double(string(char(TAB.ClippingRankSum)));
TAB.BaseQRankSum 			= double(string(char(TAB.BaseQRankSum)));
TAB.ReadPosRankSum 			= double(string(char(TAB.ReadPosRankSum)));
TAB.QD 					    = double(string(char(TAB.QD)));
TAB.MQRankSum 				= double(string(char(TAB.MQRankSum)));
TAB.PAB_MAX 				= double(string(char(TAB.PAB_MAX)));
TAB.names 					= string(TAB.names);
TAB.filter 					= string(TAB.filter);
TAB.VQSLOD 					= double(string(char(TAB.VQSLOD)));
TAB.VQSR_culprit 			= string(TAB.VQSR_culprit);
TAB.segdup 					= string(TAB.segdup);
TAB.lcr 					= string(TAB.lcr);
TAB.decoy 					= string(TAB.decoy);
TAB.nonpar 					= string(TAB.nonpar);
TAB.rf_positive_label 		= string(TAB.rf_positive_label);
TAB.rf_negative_label 		= string(TAB.rf_negative_label);
TAB.rf_label 				= string(TAB.rf_label);
TAB.rf_train 				= string(TAB.rf_train);
TAB.was_mixed 				= string(TAB.was_mixed);
TAB.has_star 				= string(TAB.has_star);
TAB.POPMAX 					= string(TAB.POPMAX);
TAB.AC_POPMAX 				= double(string(char(TAB.AC_POPMAX)));
TAB.AN_POPMAX 				= double(string(char(TAB.AN_POPMAX)));
TAB.AF_POPMAX 				= double(string(char(TAB.AF_POPMAX)));
TAB.NHOMALT_POPMAX 			= double(string(char(TAB.NHOMALT_POPMAX)));
TAB.EXON 					= string(TAB.EXON);
TAB.INTRON 					= string(TAB.INTRON);
TAB.PROTEIN_POS 			= string(TAB.PROTEIN_POS);
TAB.AMINO_ACIDS 			= string(TAB.AMINO_ACIDS);
TAB.CODONS 					= string(TAB.CODONS);
TAB.EXISTING_VAR 			= string(TAB.EXISTING_VAR);
TAB.FLAGS 					= string(TAB.FLAGS);
TAB.HGNC_ID 				= double(string(char(TAB.HGNC_ID)));
TAB.CANONICAL 				= string(TAB.CANONICAL);
TAB.CCDS 					= string(TAB.CCDS);
TAB.SWISSPROT 				= string(TAB.SWISSPROT);
TAB.TREMBL 					= string(TAB.TREMBL);
TAB.GENE_PHENO 				= string(TAB.GENE_PHENO);
TAB.HGVS_OFFSET 			= double(string(char(TAB.HGVS_OFFSET)));
TAB.CLIN_SIG 				= string(TAB.CLIN_SIG);
TAB.SOMATIC 				= string(TAB.SOMATIC);
TAB.PHENO 					= string(TAB.PHENO);
TAB.LOF_FLAGS 				= string(TAB.LOF_FLAGS);









%==========================================================================
%% EXPORT TABLE
%==========================================================================
clearvars -except P JSON TAB


save([P.gnomat P.f P.CHR '.mat'],'TAB');

