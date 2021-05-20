function [GCR] = gcr(VCR,GENEi)




% for g = 1:numel(unique(TBL.GENEi))
% 
%     TBL.GCR(TBL.GENEi==g) = 1 - prod(1 - TBL.VCR(TBL.GENEi==g));
% 
% end
% TBL = movevars(TBL,{'GCR'},'After','GENE');

GCR = nan(size(VCR));

for g = 1:numel(unique(GENEi))

    GCR(GENEi==g) = 1 - prod(1 - VCR(GENEi==g));

end


end