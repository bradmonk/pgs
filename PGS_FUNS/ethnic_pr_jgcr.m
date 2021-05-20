function [PrJGCR] = ethnic_pr_jgcr(GCR)

nB = 20;


% SQUARE AF TO GET THE JOINT PROBABILITY
GF = GCR.^2;


nV = numel(GF);
if nV > 5000
nV = 5000;
end


PrJGCR = zeros(nV,nB);


for i = 1:nV
    for j = 1:nB

        r = randperm(nV,i);

        PrJGCR(i,j) = 1 - prod(1 - GF(r));

    end
end


PrJGCR = nanmean(PrJGCR ,2);


end