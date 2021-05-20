function [VCR] = ethnicity_vcr(AN,AC,nhomalt)

VCR = (AC - nhomalt) ./ (.5 .* AN);

end