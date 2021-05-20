function [VCR] = vcr(AN,AC,nhomalt)




%TBL.VCR = (TBL.AC - TBL.nhomalt) ./ (.5 .* TBL.AN);
%TBL = movevars(TBL,{'VCR'},'After','GENE');



VCR = (AC - nhomalt) ./ (.5 .* AN);





end