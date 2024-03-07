clear all;
close all;
clc
%% 1
dt = load('../Data/SmetsWoutersData.mat');
dt.dates= datetime(dt.dates, 'ConvertFrom', 'datenum');
[T,n] = size(dt.data);

p=4;     
constant=1;
out_var=VAR_OLS(dt.data,constant,p); 

F = companionform(out_var.C);
c_tilde = [out_var.c;zeros((p-1)*n,1)];


SS = inv((eye(n)-cell2mat(out_var.C(1))-cell2mat(out_var.C(2))-cell2mat(out_var.C(3))-cell2mat(out_var.C(4))))*out_var.c;
disp("Unconditional mean VAR flat prior")
SS


hmax = 12;
start = find(dt.dates == '31-Dec-1989');

forecast_1 = zeros(hmax,n,T-start);

for j=start+1:T
    dt_sub= dt.data(1:j-1,:); %get data
  
    out_var_sub = VAR_OLS(dt_sub,true,p);
    
    F_sub = companionform(out_var_sub.C);
    c_tilde = [out_var_sub.c;zeros((p-1)*n,1)];

    yt = [dt_sub(end,:)';dt_sub(end-1,:)';dt_sub(end-2,:)';dt_sub(end-3,:)'];
    FF=F_sub;
    SF = eye(n*p);
    xihat = c_tilde + F_sub*yt;
    forecast_1(1,:,j-start)= xihat(1:n,1);
    for i=2:hmax
        SF = SF + FF;
        FF = FF*F_sub;
        xihat = SF*c_tilde + FF*yt;
        forecast_1(i,:,j-start)= xihat(1:n,1);
    end
end

%% 2. RW
forecast_2 = zeros(hmax,n,T-start);
x = 0;

for j=start+1:T
    x = x+1;
    dt_sub = dt.data(1:j-1,:); 
    for i=1:hmax
        forecast_2(i,:,x) = dt_sub(end,:); 
    end
end

%% 3. 
dt_na =  cat(1,dt.data, NaN(hmax-1,n));
dt_true = ones(hmax,n,T-start);

x = 0;
for j=start+1:T
    x = x+1;
    for i=1:hmax
        dt_true(i,:,x) = dt_na(j+i-1,:);
    end
end
 
for i=1:hmax
    rmse_1(:,:)=sqrt(nanmean((dt_true-forecast_1).^2,3));
    rmse_2(:,:)=sqrt(nanmean((dt_true-forecast_2).^2,3));
    mae_1(:,:)=nanmean(abs(dt_true-forecast_1),3);
    mae_2(:,:)=nanmean(abs(dt_true-forecast_2),3);
end

%ratios
rmse_12 = rmse_1./rmse_2;
mae_12 = mae_1./mae_2;


horizon=linspace(1,12,12)';

f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

x  = 1;

 for ii = 1:n
    subplot(2,4,x)
    
    plot(horizon, rmse_12(:,x))
    hold on;
    plot(horizon, mae_12(:,x))
    hold off;
    ylim([0.5,1.7])
    legend('RMSE', 'MAE');

    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 8); 
    set(gca,'XTick',[0:round(hmax/4.33):hmax]')
    set(gca,'XTickLabel',num2str([0:round(hmax/4.33):hmax]'))    
    title(strcat(dt.varNames{ii}))
    set(gca,'Layer','top')
    x  = x+1;
    
 end
    

saveas(gcf, '../Output/Plot/ps4_fig1.png');


%% 4. Estimate a VAR with p = 4 lags, using the entire sample
p=4; 
hmax = 12;
[T,n] = size(dt.data);

fset.prior_family = 'conjugate';
fset.prior = 'Minnesota';
fset.NoCointegration = 1;
fset.SingleUnitRoot = 1;

[~,~,out_bvar] = BVAR_v2c(dt.data,ones(T,1),p,fset);

beta = reshape(out_bvar.b_post,29,7);
[c,A,C] =reshapeVAR(beta,p,1);
    
F = companionform(C);
c_tilde = [c;zeros((p-1)*n,1)];


SS_bvar = inv((eye(n)-cell2mat(C(1))-cell2mat(C(2))-cell2mat(C(3))-cell2mat(C(4))))*c;
disp("Unconditional mean Bayesian VAR with default priors")
SS_bvar

start = find(dt.dates == '31-Dec-1989');
    
forecast_1 = zeros(hmax,n,T-start);
    
for j=start+1:T
    dt_sub = dt.data(1:j-1,:); 
    [t,n]=size(dt_sub);

    [~,~,out_bvar_sub] = BVAR_v2c(dt_sub,ones(t,1),p,fset);

    beta_sub = reshape(out_bvar_sub.b_post,29,7);
    
    [c,A,C] = reshapeVAR(beta_sub,p,1);
     F_sub = companionform(C);
     c_tilde = [c;zeros((p-1)*n,1)];
    
    yt = [dt_sub(end,:)';dt_sub(end-1,:)';dt_sub(end-2,:)';dt_sub(end-3,:)'];
    FF=F_sub;
    SF = eye(n*p);
    xihat = c_tilde + F_sub*yt;
    forecast_1(1,:,j-start)= xihat(1:n,1);
    for i=2:hmax
        SF = SF + FF;
        FF = FF*F_sub;
        xihat = SF*c_tilde + FF*yt;
        forecast_1(i,:,j-start)= xihat(1:n,1);
    end
end
    
for i=1:hmax
    rmse_bvar_1(:,:)=sqrt(nanmean((dt_true-forecast_1).^2,3));
    rmse_bvar_2(:,:)=sqrt(nanmean((dt_true-forecast_2).^2,3));
    mae_bvar_1(:,:)=nanmean(abs(dt_true-forecast_1),3);
    mae_bvar_2(:,:)=nanmean(abs(dt_true-forecast_2),3);
end


rmse_bvar_12 = rmse_bvar_1./rmse_bvar_2;
mae_bvar_12 = mae_bvar_1./mae_bvar_2;


%tighter hyperparameters
fset.alpha=0.2;
fset.lambda=0.02;
fset.mu=0.02;
fset.theta=0.1;

forecast_1 = zeros(hmax,n,T-start);
    
for j=start+1:T
    dt_sub = dt.data(1:j-1,:);
    [t,n]=size(dt_sub);

    [~,~,out_bvar_sub] = BVAR_v2c(dt_sub,ones(t,1),p,fset);

    beta_sub = reshape(out_bvar_sub.b_post,29,7);
    
    %Compute the companion form matrices 
    [c,A,C] = reshapeVAR(beta_sub,p,1);
     F_sub = companionform(C);
     c_tilde = [c;zeros((p-1)*n,1)];
    
    yt = [dt_sub(end,:)';dt_sub(end-1,:)';dt_sub(end-2,:)';dt_sub(end-3,:)'];
    FF=F_sub;
    SF = eye(n*p);
    xihat = c_tilde + F_sub*yt;
    forecast_1(1,:,j-start)= xihat(1:n,1);
    for i=2:hmax
        SF = SF + FF;
        FF = FF*F_sub;
        xihat = SF*c_tilde + FF*yt;
        forecast_1(i,:,j-start)= xihat(1:n,1);
    end
end
    
for i=1:hmax
    rmse_bvar_tight_1(:,:)=sqrt(nanmean((dt_true-forecast_1).^2,3));
    rmse_bvar_tight_2(:,:)=sqrt(nanmean((dt_true-forecast_2).^2,3));
    mae_bvar_tight_1(:,:)=nanmean(abs(dt_true-forecast_1),3);
    mae_bvar_tight_2(:,:)=nanmean(abs(dt_true-forecast_2),3);
end


rmse_bvar_tight_12 = rmse_bvar_tight_1./rmse_bvar_tight_2;
mae_bvar_tight_12 = mae_bvar_tight_1./mae_bvar_tight_2;


f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

x  = 1;
 for ii = 1:n
    subplot(2,4,x)
    
    plot(horizon, rmse_bvar_12(:,x))
    hold on;
    plot(horizon, mae_bvar_12(:,x))
    plot(horizon, rmse_bvar_tight_12(:,x))
    plot(horizon, mae_bvar_tight_12(:,x))
    hold off;
    ylim([0,1.5])
   legend('RMSE - default priors', 'MAE - default priors','RMSE - tight priors', 'MAE - tight priors');
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 8); 
    set(gca,'XTick',[0:round(hmax/4.33):hmax]')
    set(gca,'XTickLabel',num2str([0:round(hmax/4.33):hmax]'))    
    title(strcat(dt.varNames{ii}))
    set(gca,'Layer','top')
    x  = x+1;
    
 end
    

saveas(gcf, '../Output/Plot/ps4_fig2.png');

