function [PrVCR] = ethnic_pr_vcr(VCR)


nB = 20; 


nV = numel(VCR);
if nV > 5000
nV = 5000;
end


PrVCR = zeros(nV,nB);


for i = 1:nV
    for j = 1:nB

        r = randperm(nV,i);

        PrVCR(i,j) = 1 - prod(1 - VCR(r));

    end
end


PrVCR = nanmean(PrVCR,2);


end