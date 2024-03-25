function [B_draw,SIGMA_draw,posteriors] = BVAR_v2(Yraw,exogenous,p,prior_settings)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Y_{t} = c + B_{1}*Y_{t-1} + ... + B_{p}*Y_{t-p} + u_{t}

%% --------------------------DATA HANDLING---------------------------------

% Get dimensions of dependent variable matrix
[Traw, n] = size(Yraw);

% Generate lagged Y matrix. This will be part of the X matrix
Ylags = mlag2(Yraw,p); % Y is [T x M]. ylag is [T x (Mp)]

if isempty(exogenous)
    X = Ylags(p+1:Traw,:);
else
    X = [exogenous(p+1:Traw,:) Ylags(p+1:Traw,:)];
end

if size(exogenous,2)==0
    constant=0;
else constant = 1;
end

Y = Yraw(p+1:Traw,:);       % Delete first "p" rows accordingly

Y0 = mean(Yraw(1:p,:),1);   % Save the average of the first "p" rows

[T, K] = size(X);           % Get size of final matrix X


%% ------------------------- RETRIEVE PRIOR SETTINGS ----------------------

hyperparameters = [];

v2struct(prior_settings);
if exist('hyperparameters', 'var'), v2struct(hyperparameters); end
if ~exist('stationary','var'), stationary = []; end
if ~exist('iid','var'), iid = []; end

%% -------------------------DUMMY OBSERVATIONS FOR PRIOR ------------------

v0 = n+2;

Td = 0;

if exist('Y0bar','var')
    Y0(~isnan(Y0bar)) = Y0bar(~isnan(Y0bar));
end

if strcmp(prior,'Minnesota')
        
        % If parameters are not inputed use default values
         if ~exist('alpha','var'), alpha = 2; end       % lag-decaying parameter of the MN prior
         if ~exist('lambda','var'), lambda = .2; end     % std of MN prior
         if ~exist('Vc','var'), Vc = 1E6; end

         if ~exist('psi','var')                        % residual variance of AR(1) for each variable
                psi = nan(n,1);
                for i=1:n
                    ar1=ols1(Yraw(2:end,i),[ones(Traw-1,1),Yraw(1:end-1,i)]);
                    psi(i)=ar1.sig2hatols;
                end
          end

        YdummyVariance =[]; 
        XdummyVariance =[];
        for jj = 1:v0
            YdummyVariance = [YdummyVariance; diag(psi).^.5]; 
            XdummyVariance = [XdummyVariance; zeros(n,(n*p+1))]; 
        end

        Y = [YdummyVariance; Y];
        X = [XdummyVariance; X];

        Td = Td+size(YdummyVariance,1);   

        lag_1_prior = ones(1,n);           % Trending variables centered at one
        lag_1_prior(stationary) = 0.7;     % Stationary variables centered instead at zero
        lag_1_prior(iid) = 0;              % Stationary variables centered instead at zero

        % Prior on the Coefficients
        YdummyMN = [diag(lag_1_prior)*(diag(psi).^.5)/lambda; zeros(n*(p-1),n)];
        XdummyMN = [zeros(n*p,1) kron(diag([1:p])^(alpha/2),(diag(psi).^.5)/lambda)]; 

        % Prior on the Constant terms
        YdummyMN = [YdummyMN; zeros(1,n)];
        XdummyMN = [XdummyMN; [sqrt(1/Vc) zeros(1,n*p)]];

        Y = [YdummyMN; Y];
        X = [XdummyMN; X];

        Td = Td+size(YdummyMN,1);     
    
    if NoCointegration == 1     % Also known as Sum of Coefficients Prior

        if ~exist('mu','var'), mu = 1; end % If parameters are not inputed use default values
        
        % CODE UP THE SUM OF COEFFICIENTS PRIOR HERE
        y_star = zeros(n,n);
        x_star = zeros(n,n,p);

        for i=1:n
            y_star(i,i)=Y0(i)*mu^(-1);
        end
        

        for j=1:p
            for i=1:n
                x_star(i,i,j) = Y0(i)*mu^(-1);
            end
        end

        x = [x_star(:,:,1)];
        
        for i=2:p
            x = [x x_star(:,:,i)];
        end
        
        x = [zeros(n,1) x];


        X = [x;X];
        Y = [y_star;Y];

        
    end
    
    if SingleUnitRoot == 1     % Also known as Dummy Initial Observation
        
        if ~exist('theta','var'), theta = 1; end % If parameters are not inputed use default values
  
        % CODE UP THE DUMMY INITIAL OBSERVATION PRIOR HERE
        y_star = zeros(1,n);
        x_star = zeros(1,n,p);
        
        for i=1:n
            y_star(1,i) = Y0(i)*theta^(-1);
        end
        
     
        for j=1:p
            for i=1:n
                x_star(1,i,j) = Y0(i)*theta^(-1);
            end
        end

        x = [x_star(:,:,1)];
        
        for i=2:p
            x = [x x_star(:,:,i)];
        end
        
        x = [theta^(-1) x];

      
        X = [x;X];
        Y = [y_star;Y];

    end
   
    T = T + Td;
    
end
%% --------------------------OLS QUANTITIES--------------------------------

B_OLS = (X'*X)\(X'*Y); % OLS matrix of coefficients

U_OLS = Y-X*B_OLS;     % OLS residuals
SSE = U_OLS'*U_OLS;    % Residual Sum of Squares

SIGMA_OLS = SSE/T;     % OLS covariance matrix of residuals
  
%% ----------------------------- POSTERIOR ------------------------------------
                
        NT = X'*X;
        BbarT = B_OLS;

        b_post = BbarT(:);

        vT = T;
        ST = SIGMA_OLS;

        S_post = ST*vT;
        v_post = vT;
                        
%% ----------------------------- DRAW -------------------------------------
        invS_post = inv(S_post);
        
        SIGMA_draw = inv2(wish(invS_post,v_post));
        
        V_B_post = kron(SIGMA_draw,inv2(NT));
        
        V_B_post = (V_B_post + V_B_post')./2;
        
        B_draw = mvnrnd(b_post,V_B_post);
        
        B_draw = reshape(B_draw,K,n);
        
        U_post = Y-X*BbarT;     % Residuals at posterior mean

        posteriors = v2struct(invS_post,v_post,NT,vT,ST,b_post,K,n,S_post,U_post);
                
end

