function [TBL] = GENEijk(TBL,varargin)


if nargin > 1

    GLIST = varargin{1};

else

    GLIST = TBL.GENE;

end





% CREATE/SET GENEi GENEj GENEk TO ZEROS
TBL.GENEi = zeros(height(TBL),1);
TBL.GENEj = zeros(height(TBL),1);
TBL.GENEk = zeros(height(TBL),1);



% MOVE GENEi GENEj GENEk TO END OF TABLE
TBL = movevars(TBL,{'GENEi','GENEj','GENEk'},...
    'After',TBL.Properties.VariableNames{end});



% GET UNIQUE
[~,GENEj,GENEi] = unique(GLIST,'stable');



% GENEi: each gene is given a unique number from 1:Ngenes. All VARS
% for a given gene are assigned the same number

TBL.GENEi = GENEi;






% GENEj: the first VAR in each gene is given a "1"
TBL.GENEj(GENEj) = 1;
TBL.GENEj = TBL.GENEj>0;



% GENEj: the first VAR in each gene is given a unique number from 1:Ngenes
GENEk = cumsum(TBL.GENEj);
TBL.GENEk = GENEk;
TBL.GENEk(~TBL.GENEj) = 0;



TBL = movevars(TBL,{'GENEi','GENEj','GENEk'},'Before',1);




end