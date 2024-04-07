function [R2,mR2,mR2_F,R2_T,top10_s,top10_mR2] = compute_mR2(Fhat,lamhat,ve2,series)
% =========================================================================
% Compute the R-squared and marginal R-squared from
% estimated factors and factor loadings.
%
%   Fhat    = estimated factors 
%   lamhat  = factor loadings 
%   ve2     = eigenvalues of covariance matrix
%   series  = variable names
%
%   R2 = R-squared for each series for each factor
%   mR2 = marginal R-squared for each series for each factor
%   mR2_F = marginal R-squared for each factor
%   R2_T = total variation explained by all factors
%   top10_s = top 10 series that load the most on each factor
%   top10_mR2 = marginal R-squared of top10_s
% =========================================================================
% N = number of series, r = number of factors
[N,r] = size(lamhat); 

% Preallocate memory for output 
R2 = NaN(N,r);                           
mR2 = NaN(N,r);
top10_s=cell(10,r);
top10_mR2=NaN(10,r);

% Compute R-squared and marginal R-squared for each series for each factor
for i = 1:r
    R2(:,i)  = (var(Fhat(:,1:i)*lamhat(:,1:i)'))';  
    mR2(:,i) = (var(Fhat(:,i)*lamhat(:,i)'))';
end

% Compute marginal R-squared for each factor 
mR2_F = ve2./sum(ve2);
mR2_F = mR2_F(1:r)';

% Compute total variation explained by all factors
R2_T=sum(mR2_F);

% Sort series by marginal R-squared for each factor
[vals,ind] = sort(mR2,'descend');

% Get top 10 series that load the most on each factor
for i=1:r
    top10_s(:,i)=series(ind(1:10,i));
    top10_mR2(:,i)=vals(1:10,i);
end


end
