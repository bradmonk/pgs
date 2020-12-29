function [PrVCR] = pr_vcr(TBL,varargin)


if nargin == 2

    nB = varargin{1};

else

    nB = 20;

end



VCR = TBL.VCR;


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





% [~,~] = fig1ax('Number of LoF sites considered',...
%                'P(cumulative VCR)');
% scatter(1:nV, nanmean(CVCR,2), 600, [.9 .1 .1],...
%     'Marker','.','MarkerEdgeAlpha',.2);


end