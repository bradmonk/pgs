function [] = VoM(AN,AC,HOM)


PEOPLE = AN / 2;

NCARRY = AC - HOM;

VCR = NCARRY / PEOPLE;

MAF = AC / AN;

VoM = VCR / MAF - 1;

pHOM = HOM*2/AC;


fprintf('\n-----------------------------------------------------\n')
fprintf('AN     = %4.0f    CHROMOSOMES SEQUENCED\n',AN)
fprintf('AC     = %4.0f    ALT ALLELES FOUND\n',AC)
fprintf('HOM    = %4.0f    HOMOZYGOUS ALT/ALT PEOPLE\n',HOM)
fprintf('------------------------------------------\n')
fprintf('PEOPLE = AN / 2 \t\t\t= %4.0f  \n', PEOPLE)
fprintf('NCARRY = AC - HOM \t\t\t= %4.0f  \n', NCARRY)
fprintf('VCR    = NCARRY / PEOPLE \t= %4.2f  \n', VCR)
fprintf('MAF    = AC / AN \t\t\t= %4.2f  \n', MAF)
fprintf('pHOM   = HOM*2/AC \t\t\t= %4.2f  \n', pHOM)
fprintf('VoM    = VCR / MAF - 1 \t\t= %4.2f  \n', VoM)

fprintf('\n')



end