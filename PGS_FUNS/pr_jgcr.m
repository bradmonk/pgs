function [PrJGCR] = pr_jgcr(TBL, varargin)

if nargin == 2

    nB = varargin{1};

else

    nB = 20;

end





GCR = TBL.GCR(TBL.GENEj == 1).^2;

%GCR = TBL.GCR(TBL.GENEj == 1);


nV = numel(GCR);
if nV > 5000
nV = 5000;
end


PrJGCR = zeros(nV,nB);




for i = 1:nV
    for j = 1:nB

        r = randperm(nV,i);

        PrJGCR(i,j) = 1 - prod(1 - GCR(r));

    end
end



end