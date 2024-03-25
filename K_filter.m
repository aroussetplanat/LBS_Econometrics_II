function [xitt,xittm,Ptt,Pttm,loglik]=K_filter(initx,initV,x,A,C,R,Q)
% INPUTS
% x(:,t) - the observation at time t
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% initx - the initial state (column) vector 
% initV - the initial state covariance 
% OUTPUT:
% xittm = E[X(:,t) | y(:,1:t-1)]
% Pttm = Cov[X(:,t) | y(:,1:t-1)]
% xitt = E[X(:,t) | y(:,1:t)]
% Ptt = Cov[X(:,t) | y(:,1:t)]
%loglik - value of the loglikelihood

[T,N]=size(x);
r=size(A,1);

y=x';


xittm=[initx zeros(r,T)];
xitt=zeros(r,T);

Pttm=zeros(r,r,T);
Pttm(:,:,1)=initV;
Ptt=zeros(r,r,T);

for j=1:T
    
    L=inv(C*Pttm(:,:,j)*C'+R);
    nv=size(Pttm,1);
    xitt(:,j)=xittm(:,j)+Pttm(:,:,j)*C'*L*(y(:,j)-C*xittm(:,j));
    Ptt(:,:,j)=Pttm(:,:,j)-Pttm(:,:,j)*C'*L*C*Pttm(:,:,j); 
    xittm(:,j+1)=A*xitt(:,j);
    Pttm(:,:,j+1)=A*Ptt(:,:,j)*A'+Q;
    lik(j)=((2*pi)^(-N/2))*(abs((det(C*Pttm(:,:,j)*C'+R)))^(-.5))*...
        exp(-1/2*(y(:,j)-C*xittm(:,j))'*L*(-1/2*(y(:,j)-C*xittm(:,j))));
    
    
    e = y(:,j) - C*xittm(:,j); % error (innovation)
    n = length(e);
    ss = length(A);
    d = size(e,1);
    S = C*Pttm(:,:,j)*C' + R;
    GG = C'*diag(1./diag(R))*C;
    Sinv = inv(S);
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    detS = prod(diag(R))*det(eye(ss)+Pttm(:,:,j)*GG);
    denom = (2*pi)^(d/2)*sqrt(abs(detS));
    mahal = sum(e'*Sinv*e,2);
    logl(j) = -0.5*mahal - log(denom);

    loglik=sum(logl);
    
end