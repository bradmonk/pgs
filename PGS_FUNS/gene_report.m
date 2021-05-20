function [] = gene_report(T)




VNAMES = string(T.Properties.VariableNames)';

hasSV = any(VNAMES == "isSV");




if hasSV


    nGENE = sum(T.GENEj);
    nSNV  = sum(T.isSV ~=1);
    nSV   = sum(T.isSV ==1);


    fprintf('\n Genes: %4.0f  |  SNVs: %5.0f  |  SVs: %5.0f \n', nGENE , nSNV, nSV);


else


    nGENE = sum(T.GENEj);
    nSNV  = height(T);

    fprintf('\n Genes: %4.0f  |  SNVs: %5.0f \n', nGENE , nSNV );

end





end