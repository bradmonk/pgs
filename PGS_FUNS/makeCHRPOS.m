function [CHRPOS] = makeCHRPOS(CHR,POS)


    CHRPOS = uint64(CHR .* 10000000000 + POS);


end







