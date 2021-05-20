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
clc; clearvars -except P PGS


PGS = load([P.mat P.f 'PGS_STEP043_OUTPUT.mat']);

GNOMAD  = PGS.GNOMAD;





%==========================================================================
%% CREATE "RC" COLUMN FOR REF COUNTS
%==========================================================================
clc; clearvars -except P PGS GNOMAD



GNOMAD.RC = GNOMAD.AN - GNOMAD.AC;




GNOMAD = movevars(GNOMAD,{'RC'},'After','AN');
GNOMAD = hint(GNOMAD,'RC','ref allele count (+1 ref/alt; +2 ref/ref)');
hint(GNOMAD,'RC');





%==========================================================================
%% REORGANIZE TABLE COLUMNS
%==========================================================================
clc; clearvars -except P PGS GNOMAD




% GNOMAD = movevars(GNOMAD,{'geneID','rs','refBase','altBase','altAllele',},'After','AF');
% GNOMAD = movevars(GNOMAD,{'AN_raw','AN_male','AN_female'},'After','AF');
% GNOMAD = movevars(GNOMAD,{'NHOMALT'},'After','AF');





%==========================================================================
%% REMOVE SITES WITH FEW SEQUENCED ALLELES (AN)
%==========================================================================
% clc; clearvars -except P PGS GNOMAD
% 
% 
% GNOM = GNOMAD;
% 
% close all; histogram(GNOM.AN)
% 
% 
% max(GNOM.AN)            % 251,496
% min(GNOM.AN)            %       8
% median(GNOM.AN)         % 248,904
% round(mean(GNOM.AN))    % 225,813
% mode(GNOM.AN)           % 251,490
% 
% 
% 
% AN_LOWER_LIMIT = 100000; % 100k
% 
% 
% CUT = GNOM.AN < AN_LOWER_LIMIT; 
% sum( CUT)
% sum(~CUT)
% 
% GNOM( CUT , :) = [];
% 
% 
% close all; histogram(GNOM.AN)
% 
% 
% GNOMAD = GNOM;






%==========================================================================
%% SAVE A TABLE WITHOUT SEX CHROMOSOMES
%==========================================================================
% clc; clearvars -except P PGS GNOMAD
% 
% 
% 
% 
% unique(GNOMAD.CHR)
% 
% SNOMAD = GNOMAD(GNOMAD.CHR ~= 23 , :)




%==========================================================================
%% COMPUTE ALLELE FREQUENCIES
%==========================================================================
clc; clearvars -except P PGS GNOMAD
% 
% At single locus with two alleles R and A with frequencies fR and fA,
% respectively, the expected genotype frequencies under random mating are:
% 
% fRR = fR^2
% fAA = fA^2
% fRA = 2*fR*fA
%
% N.B. HARDY-WEINBERG DOES NOT WORK FOR SEX CHROMOSOMES



% GIVEN...
N   = GNOMAD.AN;
A   = GNOMAD.AC;
R   = GNOMAD.RC;
AA  = GNOMAD.nhomalt;



% COMPUTED
RA          = A - (AA .* 2);
RR          = (R - RA) ./ 2;
RR_AA_RA    = RR + AA + RA;
R_A         = R + A;
PPL         = N ./ 2;
fR          = R ./ N;
fA          = A ./ N;
fRR         = RR ./ PPL;
fAA         = AA ./ PPL;
fRA         = RA ./ PPL;

% HARDY-WEINBERG ESTIMATES
HW_fRR      = fR.^2;
HW_fAA      = fA.^2;
HW_fRA      = 2 .* fR .* fA;


FREQ = table(N, R, A, AA, RR, RA, R_A , RR_AA_RA, PPL,...
             fR, fA, fRR, fAA, fRA, HW_fRR, HW_fAA, HW_fRA);


GNOMAD.N   = N;
GNOMAD.P   = PPL;
GNOMAD.R   = R;
GNOMAD.A   = A;
GNOMAD.AA  = AA;
GNOMAD.RR  = RR;
GNOMAD.RA  = RA;
GNOMAD.fR  = fR;
GNOMAD.fA  = fA;
GNOMAD.fRR = fRR;
GNOMAD.fAA = fAA;
GNOMAD.fRA = fRA;
GNOMAD.HW_fRR = HW_fRR;
GNOMAD.HW_fAA = HW_fAA;
GNOMAD.HW_fRA = HW_fRA;





%==========================================================================
%% ADD BETTER VARIABLE DESCRIPTIONS
%==========================================================================
clc; clearvars -except P PGS GNOMAD


% GNOMAD.Properties.VariableDescriptions{'AF'} = ...
%     'AF = AC/AN; alt allele frequency';


hint(GNOMAD,'AN');
hint(GNOMAD,'AC');
hint(GNOMAD,'AF');
hint(GNOMAD,'nhomalt');
hint(GNOMAD,'VCR');
hint(GNOMAD,'GCR');

GNOMAD = hint(GNOMAD,'AN','total alleles measured (2 per person).');
GNOMAD = hint(GNOMAD,'AC','alt allele count (+1 ref|alt; +2 alt|alt)');
GNOMAD = hint(GNOMAD,'AF','AC/AN; MAF, alt allele frequency');
GNOMAD = hint(GNOMAD,'nhomalt','homozygous PEOPLE count (not chr)');
GNOMAD = hint(GNOMAD,'VCR','(AC-NHOMALT)/(AN/2); variant carrier rate');
GNOMAD = hint(GNOMAD,'GCR','1-prod(1-VCR(idx)); gene carrier rate');


GNOMAD = hint(GNOMAD,'N','Alleles measured (AN)');
GNOMAD = hint(GNOMAD,'P','People measured (N/2)');
GNOMAD = hint(GNOMAD,'A','Alt count (AC)');
GNOMAD = hint(GNOMAD,'R','Ref count (N - A)');
GNOMAD = hint(GNOMAD,'AA','Alt/Alt count (NHOMALT)');
GNOMAD = hint(GNOMAD,'RR','Ref/Ref count (R-RA)/2');
GNOMAD = hint(GNOMAD,'RA','Ref/Alt count (A - AA*2)');
GNOMAD = hint(GNOMAD,'fR','Ref frequency (R/N)');
GNOMAD = hint(GNOMAD,'fA','Alt frequency (A/N)');
GNOMAD = hint(GNOMAD,'fRR','Ref/Ref frequency (RR/P)');
GNOMAD = hint(GNOMAD,'fAA','Alt/Alt frequency (AA/P)');
GNOMAD = hint(GNOMAD,'fRA','Ref/Alt frequency (RA/P)');
GNOMAD = hint(GNOMAD,'HW_fRR','Ref/Ref Hardy-Weinberg frequency (fR^2)');
GNOMAD = hint(GNOMAD,'HW_fAA','Alt/Alt Hardy-Weinberg frequency (fA^2)');
GNOMAD = hint(GNOMAD,'HW_fRA','Ref/Alt Hardy-Weinberg frequency (2*fR*fA)');

% GNOMAD = hint(GNOMAD,'xxxxx','xxxxx');
% GNOMAD = hint(GNOMAD,'xxxxx','xxxxx');
% GNOMAD = hint(GNOMAD,'xxxxx','xxxxx');
% GNOMAD = hint(GNOMAD,'xxxxx','xxxxx');
% GNOMAD = hint(GNOMAD,'xxxxx','xxxxx');
% GNOMAD = hint(GNOMAD,'xxxxx','xxxxx');



% i = find(FREQ.fRR < 0);
% FREQ_ERROR = FREQ(i,:);
% GNO = GNOMAD(i,:);
% FREQ(GNOMAD.CHR == 23 , :) = [];


%==========================================================================
%% PLOT HARDY WEINBERG FREQUENCIES
%==========================================================================
clc; clearvars -except P PGS GNOMAD


FREQ = GNOMAD( (GNOMAD.CHR ~= 23) | (GNOMAD.CHR ~= 24) , :);



%-----------------------------------
close all;
fh1 = figure('Units','pixels','Position',[30 40 1050 700],'Color','w');
ax1 = axes('Position',[.12 .14 .82 .8],'Color','none');
xlabel('fRR & 1-fAA'); ylabel('fAA & fRR');
ARGS.XP = -.08; ARGS.YP = -.09; AXFORMAT(ax1,ARGS);
fh1.Renderer = 'painters'; hold on;
%-----------------------------------

ph1 = scatter( FREQ.fRR, FREQ.fAA ,...
    800, [.1 .1 .9],'Marker','.','MarkerEdgeAlpha',.5);
hold on;
ph2 = scatter( 1 - FREQ.fAA, FREQ.fRR ,...
    800, [.9 .1 .1],'Marker','.','MarkerEdgeAlpha',.5);
hold on;
ph3 = scatter( FREQ.fRR, FREQ.fRA ,...
    500, [.0 .3 .0],'Marker','.','MarkerEdgeAlpha',.6);
hold on;
ph4 = scatter( FREQ.fAA, FREQ.fRA ,...
    800, [.1 .9 .1],'Marker','.','MarkerEdgeAlpha',.5);




%==========================================================================
%% PLOT HARDY WEINBERG FREQUENCIES: fRR
%==========================================================================
clc; clearvars -except P PGS GNOMAD


FREQ = GNOMAD( (GNOMAD.CHR ~= 23) | (GNOMAD.CHR ~= 24) , :);



%-----------------------------------
close all;
fh1 = figure('Units','pixels','Position',[30 40 1050 700],'Color','w');
ax1 = axes('Position',[.12 .14 .82 .8],'Color','none');
xlabel('Hardy-Weinberg Estimated fRR'); ylabel('Observed fRR');
ARGS.XP = -.08; ARGS.YP = -.09; AXFORMAT(ax1,ARGS)
fh1.Renderer = 'painters'; hold on;
%-----------------------------------


ph1 = scatter( FREQ.HW_fRR, FREQ.fRR ,...
    800, [.1 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);
hold on;







%==========================================================================
%% PLOT HARDY WEINBERG FREQUENCIES: fAA
%==========================================================================
clc; clearvars -except P PGS GNOMAD


FREQ = GNOMAD( (GNOMAD.CHR ~= 23) | (GNOMAD.CHR ~= 24) , :);



%-----------------------------------
close all;
fh1 = figure('Units','pixels','Position',[30 40 1050 700],'Color','w');
ax1 = axes('Position',[.12 .14 .82 .8],'Color','none');
xlabel('Hardy-Weinberg Estimated fAA'); ylabel('Observed fAA');
ARGS.XP = -.08; ARGS.YP = -.09; AXFORMAT(ax1,ARGS)
fh1.Renderer = 'painters'; hold on;
%-----------------------------------


ph1 = scatter( FREQ.HW_fAA, FREQ.fAA ,...
    800, [.1 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);
axis([0 1 0 1])







%==========================================================================
%% PLOT HARDY WEINBERG FREQUENCIES: fRA
%==========================================================================
clc; clearvars -except P PGS GNOMAD


FREQ = GNOMAD( (GNOMAD.CHR ~= 23) | (GNOMAD.CHR ~= 24) , :);



%-----------------------------------
close all;
fh1 = figure('Units','pixels','Position',[30 40 1050 700],'Color','w');
ax1 = axes('Position',[.12 .14 .82 .8],'Color','none');
xlabel('Hardy-Weinberg Estimated fRA'); ylabel('Observed fRA');
ARGS.XP = -.08; ARGS.YP = -.09; AXFORMAT(ax1,ARGS)
fh1.Renderer = 'painters'; hold on;
%-----------------------------------


ph1 = scatter( FREQ.HW_fRA, FREQ.fRA ,...
    800, [.1 .1 .1],'Marker','.','MarkerEdgeAlpha',.4);
% axis([0 1 0 1])









%==========================================================================
%% EXPORT DATASET
%==========================================================================
clc; clearvars -except P PGS GNOMAD


DATAS = PGS.DATAS;


save([P.mat P.f 'PGS_STEP046_OUTPUT.mat'],'GNOMAD','DATAS');


%writetable(GNOMAD,[P.csv P.f 'PGS_STEP052_OUTPUT.csv'])









