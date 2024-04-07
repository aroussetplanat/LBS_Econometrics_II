clear
close all
clc
%% Q1 - 1.1
%1.Load data from CSV file
dt=importdata('../Data/current.csv',',');

% Variable names
varnames=dt.textdata(1,2:end);

% Transformation numbers
tcode=dt.data(2,:);

% SW data dummy
sw_code=dt.data(1,:);
sw_indic = find(sw_code == 1);

%Raw data
rawdata=dt.data(3:end,:);

%Dates
dates_=datetime(dt.textdata(4:end,1));

%2. Tranform variables according to tcode
data_= TransformSeries(rawdata,tcode);
%data_= prepare_missing(rawdata,tcode);

%remove first two observations as NA values for second diff transformation
data=data_(3:end,:);
dates=dates_(3:end,:);

%remove outliers (see methodo under drop_outliers)
data_nooutlier=drop_outliers(data);
%data_nooutlier=remove_outliers(data);

% sw dataset
data_nooutlier_sw = data_nooutlier(:, sw_indic);

% Estimates factors (see description in the script of the function)
r=4;
[ehat_4,Fhat_4,lamhat_4,ve2_4,data_notmi_4] = dfm_em(data_nooutlier,r);
[ehat_4_sw,Fhat_4_sw,lamhat_4_sw,ve2_4_sw,data_notmi_4_sw] = dfm_em(data_nooutlier_sw,r);


%%
%Replicate figure 1 (up to 2011:Q3)
endfig1 = find(dates == '01-Sep-2011');
r=4;
[ehat_4_fig1,Fhat_4_fig1,lamhat_4_fig1,ve2_4_fig1,~] = dfm_em(data_nooutlier(1:endfig1,:),r);
[ehat_4_sw_fig1,Fhat_4_sw_fig1,lamhat_4_sw_fig1,ve2_4_sw_fig1,~] = dfm_em(data_nooutlier_sw(1:endfig1,:),r);


f = figure;
clf
figSize = [8 5];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

for i = 1:size(Fhat_4_fig1, 2)
    subplot(size(Fhat_4_fig1, 2), 1, i);
    if i==3
        plot(dates(1:endfig1,:),-Fhat_4_fig1(:, i), 'r');
    else 
        plot(dates(1:endfig1,:),Fhat_4_fig1(:, i), 'r');
    end 
    hold on;
    if i == 1 ||i==3
        plot(dates(1:endfig1,:),-Fhat_4_sw_fig1(:, i), 'b');
    else 
        plot(dates(1:endfig1,:),Fhat_4_sw_fig1(:, i), 'b');
    end 
    hold off;    
    %ylim([-4,3])
    recessionplot
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 8); 
    title(['Factor ', num2str(i), ' Estimates']);
    legend('FRED-QD', 'S&W');
end
saveas(gcf, '../Output/Plot/final_fig1.png');

%% Q1 - 1.2 Compute the R-squared and marginal R-squared 
% see description under the function's script compute_mR2
[R2_4,mR2_4,mR2_F_4,R2_T_4,top10_s_4,top10_mR2_4] = compute_mR2(Fhat_4_fig1,lamhat_4_fig1,ve2_4_fig1,varnames);
[R2_4_sw,mR2_4_sw,mR2_F_4_sw,R2_T_4_sw,top10_s_4_sw,top10_mR2_4_sw] = compute_mR2(Fhat_4_sw_fig1,lamhat_4_sw_fig1,ve2_4_sw_fig1,varnames);


%replicate fig 2
[ehat_7,Fhat_7,lamhat_7,ve2_7,data_notmi_7] = dfm_em(data_nooutlier(1:endfig1,:),7);

 for i = 1:size(Fhat_7, 2)
 
     subplot(4, 2, 2*i - 1);
     plot(dates, Fhat_7(:, 2*i - 1), 'r');
     recessionplot
     set(gca, 'FontName', 'Times New Roman');
     set(gca, 'FontSize', 8);
     title(['Factor ', num2str(2*i - 1), ' Estimates']);

    if i == 4
        break; 
    end

     subplot(4, 2, 2*i);
     plot(dates, Fhat_7(:, 2*i), 'r');
     recessionplot
     set(gca, 'FontName', 'Times New Roman');
     set(gca, 'FontSize', 8);
     title(['Factor ', num2str(2*i), ' Estimates']);
      
 end
saveas(gcf, '../Output/Plot/final_fig2.png');


%% Q1 - 2.  Replicate all of point 1, using r = 2
r=2;
[ehat_2,Fhat_2,lamhat_2,ve2_2,data_notmi_2] = dfm_em(data_nooutlier(1:endfig1,:),r);
[ehat_2_sw,Fhat_2_sw,lamhat_2_sw,ve2_2_sw,data_notmi_2_sw] = dfm_em(data_nooutlier_sw(1:endfig1,:),r);

%Replicate figure 1
f = figure;
clf
figSize = [8 5];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

for i = 1:size(Fhat_2, 2)

    subplot(size(Fhat_2, 2), 1, i);
   
    plot(dates(1:endfig1,:),Fhat_2(:, i), 'r');
    hold on;
    if i == 1
        plot(dates(1:endfig1,:),-Fhat_2_sw(:, i), 'b');
    else
        plot(dates(1:endfig1,:),Fhat_2_sw(:, i), 'b');
    end
    hold off;
    %ylim([-4,3])
    recessionplot
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 8); 
    title(['Factor ', num2str(i), ' Estimates']);
    legend('FRED-QD', 'S&W');
end
saveas(gcf, '../Output/Plot/final_fig3.png');

%R-squared and marginal R-squared with r=2
[R2_2,mR2_2,mR2_F_2,R2_T_2,top10_s_2,top10_mR2_2] = compute_mR2(Fhat_2,lamhat_2,ve2_2,varnames);
[R2_2_sw,mR2_2_sw,mR2_F_2_sw,R2_T_2_sw,top10_s_2_sw,top10_mR2_2_sw] = compute_mR2(Fhat_2_sw,lamhat_2,ve2_2,varnames);

%% Q2. Repeat Question 1. without removing the outliers and the missing-data series.
%data_: with outliers and missing data
%not sure what to do. maybe estimate DFM via kalman
%r=4;


%% Q3.

%3.1 Rolling window forecast with factors 
% Forecast with factors (rolling window of n=102)
%number of autoregressive lags p=4 for Y and last known value of the factor
id_IP = find(strcmp(varnames, 'INDPRO'));
id_CPI = find(strcmp(varnames, 'CPIAUCSL'));
data_sub = [data_nooutlier(:,id_IP),data_nooutlier(:,id_CPI)];
[T,n]=size(data_sub);
start = find(dates == '01-Mar-1985');

hmax = [1,4,8]; % forecast horizon
p=4;
window_size= start-1;
forecast_1 = NaN(length(hmax),n,T-2-start); %with factors from 1.1
forecast_2 = NaN(length(hmax),n,T-2-start); %with factors from 1.2
forecast_3 = NaN(length(hmax),n,T-2-start); %with factors from 2.1
forecast_4 = NaN(length(hmax),n,T-2-start); %with factors from 2.2
forecast_true = NaN(length(hmax),n,T-2-start);

for j=start+1:(T-2) %loop over window

    % not outliers
    data_nooutlier_window=data_nooutlier(j-window_size-1:j-1,:);
    x_na = mean(isnan(data_nooutlier_window)) > 0.99;
    data_nooutlier_window(:, x_na) = [];
            
    %keeping outliers
    data_window=data(j-window_size-1:j-1,:);
    x_na = mean(isnan(data_window)) > 0.99;
    data_window(:, x_na) = [];
     
    %re-estimate factors with available info:
    [~,Fhat_4_nooutlier_window,~,~,~] = dfm_em(data_nooutlier_window,4);
    [~,Fhat_2_nooutlier_window,~,~,~] = dfm_em(data_nooutlier_window,2);
    [~,Fhat_4_window,~,~,~] = dfm_em(data_window,4);
    [~,Fhat_2_window,~,~,~] = dfm_em(data_window,2);
    
    for v= 1:n %loop over variables
            y_window = data_sub(j-window_size-1:j-1,v); 

        for i = 1:length(hmax) %loop over horizon forecast
            h = hmax(i);

            if j+h  > T-1
                break
            end
            y_window_lag=lagmatrix(y_window, h:(h+p-1));
            Fhat_4_nooutlier_window_lag=lagmatrix(Fhat_4_nooutlier_window, h);
            Fhat_2_nooutlier_window_lag=lagmatrix(Fhat_2_nooutlier_window, h);
            Fhat_4_window_lag=lagmatrix(Fhat_4_window, h);
            Fhat_2_window_lag=lagmatrix(Fhat_2_window, h);

            %Forecast 1: 4 factors no outliers
            X = [ones(length(y_window), 1), y_window_lag, Fhat_4_nooutlier_window_lag];
            beta_hat = regress(y_window, X);
            forecast_1(i,v,j-start-1+h) =beta_hat(1)+beta_hat(2)*y_window(end)+beta_hat(3)*y_window(end-1)+beta_hat(4)*y_window(end-2)+beta_hat(5)*y_window(end-3)+Fhat_4_nooutlier_window(end,:)*beta_hat(6:end);

            %Forecast 2: 2 factors no outliers
            X = [ones(length(y_window), 1), y_window_lag, Fhat_2_nooutlier_window_lag];
            beta_hat = regress(y_window, X);
            forecast_2(i,v,j-start-1+h) =beta_hat(1)+beta_hat(2)*y_window(end)+beta_hat(3)*y_window(end-1)+beta_hat(4)*y_window(end-2)+beta_hat(5)*y_window(end-3)+Fhat_2_nooutlier_window(end,:)*beta_hat(6:end);

             %Forecast 3: 4 factors with outliers
            X = [ones(length(y_window), 1), y_window_lag, Fhat_4_window_lag];
            beta_hat = regress(y_window, X);
            forecast_3(i,v,j-start-1+h) =beta_hat(1)+beta_hat(2)*y_window(end)+beta_hat(3)*y_window(end-1)+beta_hat(4)*y_window(end-2)+beta_hat(5)*y_window(end-3)+Fhat_4_window(end,:)*beta_hat(6:end);

            %Forecast 2: 2 factors with outliers
            X = [ones(length(y_window), 1), y_window_lag, Fhat_2_window_lag];
            beta_hat = regress(y_window, X);
            forecast_4(i,v,j-start-1+h) =beta_hat(1)+beta_hat(2)*y_window(end)+beta_hat(3)*y_window(end-1)+beta_hat(4)*y_window(end-2)+beta_hat(5)*y_window(end-3)+Fhat_2_window(end,:)*beta_hat(6:end);

            %true
            forecast_true(i,v,j-start-1+h)=data_sub(j-start-1+h,v);

        end 
    end
end


      
%% 3.2 Random Walk forecast 
forecast_5 = zeros(length(hmax),n,T-start-2); %rw 

for j=start+1:(T-2)
    y_window = data_sub(j-window_size-1:j-1,:); 
    mu=nanmean(y_window);
    for i = 1:length(hmax) %loop over horizon forecast
            h = hmax(i);
             if j+h  > T-1
                break
            end
        forecast_5(i,:,j-start-1+h) = mu*h; 
    end  
end

%% 3.3

for i=1:length(hmax)
    rmse_1(:,:)=sqrt(nanmean((forecast_true-forecast_1).^2,3));
    rmse_2(:,:)=sqrt(nanmean((forecast_true-forecast_2).^2,3));
    rmse_3(:,:)=sqrt(nanmean((forecast_true-forecast_3).^2,3));
    rmse_4(:,:)=sqrt(nanmean((forecast_true-forecast_4).^2,3));
    rmse_5(:,:)=sqrt(nanmean((forecast_true-forecast_5).^2,3));

    mae_1(:,:)=nanmean(abs(forecast_true-forecast_1),3);
    mae_2(:,:)=nanmean(abs(forecast_true-forecast_2),3);
    mae_3(:,:)=nanmean(abs(forecast_true-forecast_3),3);
    mae_4(:,:)=nanmean(abs(forecast_true-forecast_4),3);
    mae_5(:,:)=nanmean(abs(forecast_true-forecast_5),3);

end

%ratios wrt to random walk 
rmse_15 = rmse_1./rmse_5;
rmse_25 = rmse_2./rmse_5;
rmse_35 = rmse_3./rmse_5;
rmse_45 = rmse_4./rmse_5;

mae_15 = mae_1./mae_5;
mae_25 = mae_2./mae_5;
mae_35 = mae_3./mae_5;
mae_45 = mae_4./mae_5;


horizon=hmax';

f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])
varnames={'IP', 'CPI'};
x  = 1;

 for ii = 1:n
    subplot(2,2,x)
    
    plot(horizon, rmse_15(:,x))
    hold on;
    plot(horizon, rmse_25(:,x))
    hold on;
    plot(horizon, rmse_35(:,x))
    hold on;
    plot(horizon, rmse_45(:,x))
    hold off;
    %ylim([0.5,1.7])
    legend('4 factors', '2 factors', '4 factors with outliers and missing', '2 factors with outliers and missing');

    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 9); 
    %set(gca,'XTick',[0:round(hmax/4.33):hmax]')
    %set(gca,'XTickLabel',num2str([0:round(hmax/4.33):hmax]'))    
    title(varnames{ii},' RMSE')
    set(gca,'Layer','top')

    subplot(2,2,x+2)
    plot(horizon, mae_15(:,x))
    hold on;
    plot(horizon, mae_25(:,x))
    hold on;
    plot(horizon, mae_35(:,x))
    hold on;
    plot(horizon, mae_45(:,x))
    hold off;
    %ylim([0.5,1.7])
    legend('4 factors', '2 factors', '4 factors with outliers and missing', '2 factors with outliers and missing');

    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 9); 
    set(gca,'XTick',[0:round(hmax/4.33):hmax]')
    set(gca,'XTickLabel',num2str([0:round(hmax/4.33):hmax]'))    
    title(varnames{ii}, ' MAE')
    set(gca,'Layer','top')
    x  = x+1;
 end
    
saveas(gcf, '../Output/Plot/final_fig4.png');

%% Q4 - 1
p=4; 
fset.prior = 'Minnesota';
fset.NoCointegration = 1;
fset.SingleUnitRoot = 1;
fset.alpha = 2;
fset.lambda = 0.2;
fset.mu = 1;
fset.theta = 1;

x_na=isnan(data_sub);
mu=repmat(nanmean(data_sub),T,1);

% Replace missing values with unconditional mean
x_notmi=data_sub;
x_notmi(x_na)=mu(x_na);

forecast_6 = NaN(length(hmax),2,T-start-2); %4 factors
forecast_7 = NaN(length(hmax),2,T-start-2); %2 factors

for j = start+1:(T-2) 
    
    y_window = x_notmi(j-window_size-1:j-1, :); 

    data_nooutlier_window=data_nooutlier(j-window_size-1:j-1,:);
    x_na = mean(isnan(data_nooutlier_window)) > 0.99;
    data_nooutlier_window(:, x_na) = [];
             
    %re-estimate factors with available info:
    [~,Fhat_4_nooutlier_window,~,~,~] = dfm_em(data_nooutlier_window,4);
    [~,Fhat_2_nooutlier_window,~,~,~] = dfm_em(data_nooutlier_window,2);
     
    % BVAR 4 factors
    Y = [y_window, Fhat_4_nooutlier_window];
    [TT, n] = size(Y);

    [~,~,out_bvar_sub] = BVAR_v2c(Y, ones(TT,1), p, fset);
    Beta = reshape(out_bvar_sub.b_post, 1 + n*p, n); 
    A = Beta(2:end,:)'; 
    c0 = Beta(1, :)';

    [n, np] = size(A);
    p = np / n;

    idx = [1:n:n*p, n*p+1];
    for u = 1:p
        C{u} = A(1:n, idx(u):(idx(u+1)-1));
    end

    % Forecast
    y_sample = Y(end-3:end, :); % four lagged periods to predict
    for h = 1:max(hmax)
        if j+h  > T-1
            break
        end
        yhat = c0 + C{1} * y_sample(4, :)' + C{2} * y_sample(3, :)' + C{3} * y_sample(2, :)' + C{4} * y_sample(1, :)';
        idx = find(hmax == h);

        if ~isempty(idx)
            forecast_6(idx,1,j-start-1+h) = yhat(1);
            forecast_6(idx,2,j-start-1+h) = yhat(2);
        end

        y_sample = [y_sample; yhat'];
        y_sample = y_sample(2:5, :);
    end

  
    % BVAR 2 factors
    Y = [y_window, Fhat_2_nooutlier_window];
    [TT, n] = size(Y);

     [~,~,out_bvar_sub] = BVAR_v2c(Y, ones(TT,1), p, fset);
    Beta = reshape(out_bvar_sub.b_post, 1 + n*p, n); 
    A = Beta(2:end,:)'; 
    c0 = Beta(1, :)';

    [n, np] = size(A);
    p = np / n;

    idx = [1:n:n*p, n*p+1];
    for u = 1:p
        C{u} = A(1:n, idx(u):(idx(u+1)-1));
    end

    % Forecast
    y_sample = Y(end-3:end, :); % four lagged periods to predict
    for h = 1:max(hmax)
        if j+h  > T-1
            break
        end
        yhat = c0 + C{1} * y_sample(4, :)' + C{2} * y_sample(3, :)' + C{3} * y_sample(2, :)' + C{4} * y_sample(1, :)';
        idx = find(hmax == h);

        if ~isempty(idx)
            forecast_7(idx,1,j-start-1+h) = yhat(1);
            forecast_7(idx,2,j-start-1+h) = yhat(2);
        end

        y_sample = [y_sample; yhat'];
        y_sample = y_sample(2:5, :);
    end

end
    

rmse_6(:,:)=sqrt(nanmean((forecast_true-forecast_6(:,1:2,:)).^2,3));
mae_6(:,:)=nanmean(abs(forecast_true-forecast_6(:,1:2,:)),3);

rmse_7(:,:)=sqrt(nanmean((forecast_true-forecast_7(:,1:2,:)).^2,3));
mae_7(:,:)=nanmean(abs(forecast_true-forecast_7(:,1:2,:)),3);

rmse_65 = rmse_6./rmse_5;
mae_65 = mae_6./mae_5;

rmse_75 = rmse_7./rmse_5;
mae_75 = mae_7./mae_5;

f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

horizon=[1 4 8];
x=1;
 for ii = 1:2

    subplot(2,2,x)
    plot(horizon, rmse_15(:,x))
    hold on;
    plot(horizon, rmse_25(:,x))
    hold on;
    plot(horizon, rmse_35(:,x))
    hold on;
    plot(horizon, rmse_45(:,x))
    hold on;
    plot(horizon, rmse_65(:,x))
    hold on;
    plot(horizon, rmse_75(:,x))
    hold off;
    %ylim([0.5,1.7])
    legend('4 factors', '2 factors', '4 factors with outliers and missing', '2 factors with outliers and missing', 'BVAR - 4 factors', 'BVAR - 2 factors');

    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 9); 
    %set(gca,'XTick',[0:round(hmax/4.33):hmax]')
    %set(gca,'XTickLabel',num2str([0:round(hmax/4.33):hmax]'))    
    title(varnames{ii},' RMSE')
    set(gca,'Layer','top')

    subplot(2,2,x+2)
    plot(horizon, mae_15(:,x))
    hold on;
    plot(horizon, mae_25(:,x))
    hold on;
    plot(horizon, mae_35(:,x))
    hold on;
    plot(horizon, mae_45(:,x))
    hold on;
    plot(horizon, mae_65(:,x))
    hold on;
    plot(horizon, mae_75(:,x))
    hold off;
    %ylim([0.5,1.7])
    legend('4 factors', '2 factors', '4 factors with outliers and missing', '2 factors with outliers and missing', 'BVAR - 4 factors', 'BVAR - 2 factors');

    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 9); 
    set(gca,'XTick',[0:round(hmax/4.33):hmax]')
    set(gca,'XTickLabel',num2str([0:round(hmax/4.33):hmax]'))    
    title(varnames{ii}, ' MAE')
    set(gca,'Layer','top')
    x  = x+1;
 end
    
saveas(gcf, '../Output/Plot/final_fig5.png');


