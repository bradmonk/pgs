function [PrJVCR] = pr_jvcr(TBL,varargin)

if nargin == 2

    nB = varargin{1};

else

    nB = 20;

end



VCR = TBL.VCR;


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





end