% getFEVD.m - Obtaining Forecast Error Variance Decomposition
function FEVD = getFEVD(IRF,hmax)
    
     n = size(IRF, 1);
    denom = sum(IRF.^2, 2);
    FEVD = zeros(n,n,hmax+1);
    
    for var = 1:n
        for h = 0:hmax
            FEVD(var,:, h+1) = sum(IRF(var,:,1:h+1).^2, 3) ./ sum(denom(var,:,1:h+1), 3);
        end
    
    end 
   
end
