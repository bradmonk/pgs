function [TBL] = vcr(TBL)

TBL.VCR = (TBL.AC - TBL.NHOMALT) ./ (.5 .* TBL.AN);


TBL = movevars(TBL,{'VCR'},'After','GENE');

end