function [PrGCR] = ethnic_pr_gcr(GCR)
%GCR = TBL.GCR(TBL.GENEj == 1);

nB = 20;


nV = numel(GCR);
if nV > 5000
nV = 5000;
end


PrGCR = zeros(nV,nB);


for i = 1:nV
    for j = 1:nB

        r = randperm(nV,i);

        PrGCR(i,j) = 1 - prod(1 - GCR(r));

    end
end


PrGCR = nanmean(PrGCR,2);


end