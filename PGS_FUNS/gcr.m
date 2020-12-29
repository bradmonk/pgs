function [TBL] = gcr(TBL)




for g = 1:numel(unique(TBL.GENEi))

    TBL.GCR(TBL.GENEi==g) = 1 - prod(1 - TBL.VCR(TBL.GENEi==g));

end


TBL = movevars(TBL,{'GCR'},'After','GENE');



end