clear all;
close all;
clc
%% 1.1.2
rng(1);

T=500;
mu=1;
mu_t = 0;
phi=0.7;
theta=0.3;
sigma=1;
y = zeros(T,1);
eps = normrnd(0,sigma,T,1);

beta = [1 theta];
R = 0;
F = [phi, 1 ; 0, 0];
H = [1,0];
Q = [sigma^2, theta*sigma; theta*sigma, theta^2*sigma^2];
A = mu-mu_t;


ep_t = [eps, theta*eps]';
y_ss=zeros(T,1);
x = [zeros(T,1) theta*eps]';
y_ss(1) = mu;
x(1,1)=y_ss(1,1)-mu;

for i=2:T
    x(:,i) = F*x(:,i-1)+ ep_t(:,i);
    y_ss(i) = A + H*x(:,i) + mu_t;
end

f = figure;
clf
figSize = [8 4];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
plot(1:T, y_ss);
saveas(gcf,'../Output/Plot/ps3_fig1.png')

%% 1.3.5
a1 = [0 0]';
p1 = (1/(1-phi^2))*eye(2);
[x_t, x_t_1, sigma_t, sigma_t_1, L] = KalmanFilter(y_ss, F,Q,A,H,R,a1,p1);

%% 1.4
est_gamma = fminunc(@(gamma) kf_lik(gamma,y_ss,a1,p1,@KalmanFilter), [mu,phi,theta,sigma]);
[L_est, x_t_est] = kf_lik(est_gamma,y_ss,a1,p1,@KalmanFilter);

%% 1.5
f = figure;
clf
figSize = [8 4];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
plot(1:T,y_ss, 1:T, x(1,:)',1:T, x_t_est(:,1));
legend({'Y_t','x_t', 'x_{t|t}'});
saveas(gcf,'../Output/Plot/ps3_fig2.png');

%% 2.1
for i=200:250
    y_ss(i,:)=NaN;
end

for i=412:460
    y_ss(i,:)=NaN;
end

%% 2.2
est_gamma_mi = fminunc(@(gam) kf_lik(gam,y_ss,a1,p1,@KalmanFilter_missingdata), [mu,phi,theta,sigma]);
[L_est_mi, x_t_est_mi] = kf_lik(est_gamma_mi,y_ss,a1,p1,@KalmanFilter_missingdata);

f = figure;
clf
figSize = [8 4];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
plot(1:T,y_ss, 1:T, x(1,:)',1:T, x_t_est_mi(:,1));
legend({'Y_t','x_t', 'x_{t|t}'});
saveas(gcf,'../Output/Plot/ps3_fig3.png');


