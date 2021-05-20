function [PrJVCR] = ethnic_pr_jvcr(VCR)

nB = 20;


% SQUARE AF TO GET THE JOINT PROBABILITY
AF = VCR.^2;


nV = numel(AF);
if nV > 5000
nV = 5000;
end


PrJVCR = zeros(nV,nB);


for i = 1:nV
    for j = 1:nB

        r = randperm(nV,i);

        PrJVCR(i,j) = 1 - prod(1 - AF(r));

    end
end


PrJVCR = nanmean(PrJVCR ,2);


end