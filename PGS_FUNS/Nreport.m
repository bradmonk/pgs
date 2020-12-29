function [] = Nreport(TBL1,varargin)


if nargin == 1

    fprintf('\n\nTABLE ROWS: %7.0f \n\n', height(TBL1))

elseif nargin == 2
    TBL2 = varargin{1};
    fprintf('\n\nTABLE-1 ROWS: %7.0f\nTABLE-2 ROWS: %7.0f\n\n', height(TBL1), height(TBL2))

elseif nargin == 3
    TBL2 = varargin{1};
    TBL3 = varargin{2};
    fprintf('\n\nTABLE-1 ROWS: %7.0f\nTABLE-2 ROWS: %7.0f\nTABLE-3 ROWS: %7.0f\n\n', height(TBL1), height(TBL2), height(TBL3))

elseif nargin == 4
    TBL2 = varargin{1};
    TBL3 = varargin{2};
    TBL4 = varargin{3};
    fprintf('\n\nTABLE-1 ROWS: %7.0f\nTABLE-2 ROWS: %7.0f\nTABLE-3 ROWS: %7.0f\nTABLE-4 ROWS: %7.0f\n\n',...
    height(TBL1), height(TBL2), height(TBL3), height(TBL4))

elseif nargin == 5
    TBL2 = varargin{1};
    TBL3 = varargin{2};
    TBL4 = varargin{3};
    TBL5 = varargin{4};
    fprintf('\n\nTABLE-1 ROWS: %7.0f\nTABLE-2 ROWS: %7.0f\nTABLE-3 ROWS: %7.0f\nTABLE-4 ROWS: %7.0f\nTABLE-5 ROWS: %7.0f\n\n',...
    height(TBL1), height(TBL2), height(TBL3), height(TBL4), height(TBL5))
end






end