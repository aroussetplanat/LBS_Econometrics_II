function [ehat,Fhat,lamhat,ve2,x_notmi]=dfm_em(x,r)
% -------------------------------------------------------------------------
%Parameters 

% Maximum number of iterations for the EM algorithm
maxit=200;

% Number of observations per series and number of series
[T,N]=size(x);

% Set error
err=2000;

% Iteration counter
it=0;

% -------------------------------------------------------------------------
% EM ALGORITHM -- Initialisation
% 1. Fill in missing values for each series with the unconditional mean of
% that series if x_na==1. Otherwise, do not replace missing values.
% 2. Demean and standardize the dataset. 
% 3. Estimate factors and use them to predict the original dataset.

%get missing values
x_na=isnan(x);

% Get unconditional mean of the non-missing values of each series
mu=repmat(nanmean(x),T,1);
% Replace missing values with unconditional mean
x_notmi=x;
x_notmi(x_na)=mu(x_na);

% Demean and standardise
mu=mean(x_notmi);
sd=std(x_notmi);
x_std=(x_notmi-mu)./sd;

% Compute principal components on std dataset using pc() (see description
% at the end of the script)
%   xhat   = values of the original dataset (i.e., x_notmi) predicted by the factors
%   Fhat   = factors scaled by (1/sqrt(N)) (N: number of series)
%   lamhat = factor loadings scaled by N
%   ve2    = eigenvalues of x_std'*x3 

[xhat,Fhat,lamhat,ve2]  = pc(x_std,r);

% Save predicted series values
xhat0=xhat;

% -------------------------------------------------------------------------
% EM ALGORITHM -- Iteration until convergence
% 1. Replace missing observations using the predicted series values from the latest estimated factors 
% 2. Demean and standardize the dataset
% 3. Estimate the factors
% 4. Repeat the process until the factor estimates do not change (i.e. diff
% between old and new factors estimated is less than err=0.0000001) or when
% maxit=50 is reached

while err> 0.000001 && it <maxit
    it=it+1;
    
    fprintf('Iteration %d: obj %10f r: %d \n',it,err,r);

    % Replace missing obs with latest predictions from the factors
    for t=1:T
        for j=1:N
            if x_na(t,j)==1
                x_notmi(t,j)=xhat(t,j)*sd(j)+mu(j);    
            else
                x_notmi(t,j)=x(t,j);
            end
        end
    end
    
    %Standardise
    mu=mean(x_notmi);
    sd=std(x_notmi);
    x_std=(x_notmi-mu)./sd;

    % Run principal components on the new dataset using pc()
    [xhat,Fhat,lamhat,ve2]  = pc(x_std,r);

    % Compute diff between the predicted values of the new dataset
    % and the predicted values of the previous dataset
    diff=xhat-xhat0;
    
    % The error value: sum of the squared differences
    % between xhat and xhat0 divided by sum of the xhat0 squared
    a1=diff(:);
    a2=xhat0(:);
    err=(a1'*a1)/(a2'*a2);

    % Set xhat0 to latest xhat
    xhat0=xhat;
end

% Calculate the diff between the initial dataset and the predicted values
% from the final factors
ehat = x-xhat.*sd-mu;


 function [xhat,fhat,lambda,ss]=pc(X,nfac)
% -------------------------------------------------------------------------
% Run principal component analysis.
% -------------------------------------------------------------------------
%           xhat   = values of X predicted by the factors
%           fhat   = factors scaled by (1/sqrt(N))  (N: number of variables)
%           lambda = factor loadings scaled by sqrt(N)
%           ss     = eigenvalues of X'*X 
% -------------------------------------------------------------------------

% Number of series
N=size(X,2);

% Singular value decomposition: X'*X = U*S*V'
[U,S,~]=svd(X'*X);

% Factor loadings scaled by sqrt(N)
lambda=U(:,1:nfac)*sqrt(N);

% Factors scaled by 1/(N) 
fhat=X*lambda/N;

% Estimate initial dataset X using the factors
xhat=fhat*lambda';

% Eigenvalues of X'*X
ss=diag(S);

end



end

