%% 1.1 A Monte Carlo exercise 
%Q1
phi = 0.9;
T = 280;
M = 5000;
unconditional_mean = 3;
unconditional_std = 3;
rng(42); 

mu = unconditional_mean * (1 - phi);
eps_sigma2 = (unconditional_std^2) * (1 - phi^2);
y = zeros(M, T);

for i = 1:M
    epsilon = sqrt(eps_sigma2) * randn(1, T);
    y(i, 1) = mu/(1 - phi) + epsilon(1);
    
    for t = 2:T
        y(i, t) = mu + phi * y(i, t-1) + epsilon(t);
    end
end


for i = 1:M
    color = rand(1, 3); 
    plot(1:T, y(i, :), 'Color', color, 'LineWidth', 1);
    hold on;
end


mean_y = mean(y);
std_y = std(y);
mean_line = plot(1:T, mean_y, 'k-', 'LineWidth', 2); 
std_plus_line = plot(1:T, mean_y + std_y, 'r--', 'LineWidth', 1.5); 
std_minus_line = plot(1:T, mean_y - std_y, 'r--', 'LineWidth', 1.5); 
%title('Simulated AR(1) Processes with Mean=3 and Standard Deviation=3');
xlabel('Time (t)');
ylabel('y_t');
grid on;
legend([mean_line, std_plus_line(1), std_minus_line(1)], 'Mean', 'Mean + Std', 'Mean - Std');

hold off;
saveas(gcf, '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Plot/ps1_fig1', 'epsc');

%% Q2

mu_hat = zeros(M, 1);
phi_hat = zeros(M, 1);

for i = 1:M
    X = [ones(T-1, 1) y(i, 1:T-1)'];
    y_t = y(i, 2:T)';
    beta_hat = X\y_t;
    mu_hat(i) = beta_hat(1);
    phi_hat(i) = beta_hat(2);
 
end


figure;
subplot(2, 1, 1);
histogram(mu_hat, 'Normalization', 'probability', 'EdgeColor', 'w');
hold on;
line([mu mu], [0 0.25], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
title('\mu OLS Estimates Distribution');
xlabel('\mu');
legend('OLS Estimates', 'True Value');

subplot(2, 1, 2);
histogram(phi_hat, 'Normalization', 'probability', 'EdgeColor', 'w');
hold on;
line([phi phi], [0 0.25], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
title('\phi OLS Estimates Distribution');
xlabel('\phi');
legend('OLS Estimates', 'True Value');

saveas(gcf, '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Plot/ps1_fig2', 'epsc');
%% Q3

M = 5000;   
phi_values = [0.90, 0.95, 0.97, 0.99];  
T_values = [40, 80, 120, 280];     
rng(42);  

unconditional_mean = 3;
unconditional_std = 3;

results_table = table();

for phi_idx = 1:length(phi_values)
    phi = phi_values(phi_idx);
    
    for T_idx = 1:length(T_values)
        T = T_values(T_idx);
      
        mu = unconditional_mean * (1 - phi);
        eps_sigma2 = (unconditional_std^2) * (1 - phi^2);

        y = zeros(M, T);

        for i = 1:M
            epsilon = sqrt(eps_sigma2) * randn(1, T);
            y(i, 1) = mu/(1 - phi) + epsilon(1);
    
            for t = 2:T
                y(i, t) = mu + phi * y(i, t-1) + epsilon(t);
            end
        end

        mu_hat = zeros(M, 1);
        phi_hat = zeros(M, 1);

    for i = 1:M
        X = [ones(T-1, 1) y(i, 1:T-1)'];
        y_t = y(i, 2:T)';
        beta_hat = X\y_t;
        mu_hat(i) = beta_hat(1);
        phi_hat(i) = beta_hat(2);
    end
        
        mean_phi_hat = mean(phi_hat);
        
        result_entry = table(phi, T, mean_phi_hat);
        results_table = [results_table; result_entry];
    end
end

disp(results_table);

latex_file = '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Table/ps1_table1.tex';
fid = fopen(latex_file, 'w');

fprintf(fid, '\\begin{table}\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{ccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'True $\\phi$ & Sample Size $T$ & Mean OLS Estimate of $\\phi$ \\\\ \n');
fprintf(fid, '\\hline\n');

for row_idx = 1:height(results_table)
    fprintf(fid, '%.2f & %d & %.4f \\\\ \n', results_table.phi(row_idx), results_table.T(row_idx), results_table.mean_phi_hat(row_idx));
end


fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Mean OLS Estimate of $\\phi$ for Different True $\\phi$ Values and Sample Sizes}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);

figure;
colors = lines(length(phi_values));
for phi_idx = 1:length(phi_values)
    phi = phi_values(phi_idx);
   
    filtered_results = results_table(results_table.phi == phi, :); 
    line_color = colors(phi_idx, :);
    subplot(2, 1, 1);  
    hold on;
    plot(filtered_results.T, filtered_results.mean_phi_hat, 'o-', 'LineWidth', 1.5, 'DisplayName', ['Estimated \phi=', num2str(phi)], 'Color', line_color);
    plot([min(filtered_results.T), max(filtered_results.T)], [phi, phi], '--', 'LineWidth', 1.5, 'Color', line_color, 'DisplayName', ['True \phi=', num2str(phi)]);
    hold off;
    diff_values = filtered_results.mean_phi_hat - phi;
    subplot(2, 1, 2);  
    hold on;
    plot(filtered_results.T, diff_values, 'o-', 'LineWidth', 1.5, 'DisplayName', ['True - Estimated \phi=', num2str(phi)], 'Color', line_color);
    hold off;
end

subplot(2, 1, 1);
xlabel('Sample Size (T)');
ylabel('Mean OLS Estimate of \phi');
title('Mean OLS Estimate of \phi for Different True \phi Values');
legend('show');

subplot(2, 1, 2);
xlabel('Sample Size (T)');
ylabel('Difference between True and Estimated \phi');
title('Difference between True and Estimated \phi for Different True \phi Values');
legend('show');
saveas(gcf, '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Plot/ps1_fig3', 'epsc');
%% 1.2 Forecasting Under Structural Breaks
%Q1
rng(8);  
phi = 0.9;
unconditional_std = 3;
eps_sigma2 = (unconditional_std^2) * (1 - phi^2);

start_year = 1948;
end_year = 2017.75;
quarters = (start_year:0.25:end_year)';

mu = zeros(size(quarters));
mu(quarters >= 1948 & quarters <= 1972) = 5 * (1 - phi);
mu(quarters > 1972 & quarters <= 1995) = 0 * (1 - phi);
mu(quarters > 1995 & quarters <= 2004) = 4.5 * (1 - phi);
mu(quarters > 2004) = 0.5 * (1 - phi);

y = zeros(size(quarters));
epsilon = sqrt(eps_sigma2) * randn(size(quarters));
y(1) = mu(1)/(1 - phi) + epsilon(1);

for t = 2:length(quarters)
    y(t) = mu(t) + phi * y(t - 1) + epsilon(t);
end


%% Q2
mean_y_1948_1972 = mean(y(quarters >= 1948 & quarters <= 1972));
mean_y_1972_1995 = mean(y(quarters > 1972 & quarters <= 1995));
mean_y_1995_2004 = mean(y(quarters > 1995 & quarters <= 2004));
mean_y_2004_2017 = mean(y(quarters > 2004));

figure;
plot(quarters, y, 'LineWidth', 1.5, 'DisplayName', 'y_t');
hold on;
plot([1948 1972], [mean_y_1948_1972 mean_y_1948_1972], 'g--', 'LineWidth', 1.5, 'DisplayName', 'Mean (1948-1972)');
plot([1972 1995], [mean_y_1972_1995 mean_y_1972_1995], 'b--', 'LineWidth', 1.5, 'DisplayName', 'Mean (1972-1995)');
plot([1995 2004], [mean_y_1995_2004 mean_y_1995_2004], 'm--', 'LineWidth', 1.5, 'DisplayName', 'Mean (1995-2004)');
plot([2004 2017], [mean_y_2004_2017 mean_y_2004_2017], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean (2004-2017)');

xlabel('Year');
legend('show');
saveas(gcf, '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Plot/ps1_fig4', 'epsc');

ytrue = y;
%% Q3
start_forecast = find(quarters == 1990.00);

forecast_matrix_1 = zeros(12, length(quarters) - start_forecast+1);

for i = start_forecast:length(quarters)
    y_window = ytrue(1:i-2);

    X = [ones(length(y_window), 1), y_window];
    y = ytrue(2:i-1); 
    beta_hat = X\y;

    forecast_horizon = 12;
    forecasts = zeros(forecast_horizon, 1);
    for h = 1:forecast_horizon
        forecasts(h) = beta_hat(1) + beta_hat(2)^(h) * (y(end) - beta_hat(1)) ;
    end

     forecast_matrix_1(:, i-start_forecast+1) = forecasts;
end

%% Q4
forecast_matrix_2 = zeros(12, length(quarters) - start_forecast+1);
rolling_window_size = 40;

for i = start_forecast:length(quarters)
    yt_window = ytrue(i - rolling_window_size - 2:i-2);

    X = [ones(length(yt_window), 1), yt_window];
    y = ytrue(i - rolling_window_size-1:i-1); 
    beta_hat = X\y;

    forecast_horizon = 12;
    forecasts = zeros(forecast_horizon, 1);
    for h = 1:forecast_horizon
        forecasts(h) = beta_hat(1) + beta_hat(2)^(h) * (y(end) - beta_hat(1));
    end

    forecast_matrix_2(:, i-start_forecast+1) = forecasts;
end

%% Q5
forecast_matrix_3 = zeros(12, length(quarters) - start_forecast+1);


for i = start_forecast:length(quarters)
    yt_window = ytrue(1:i-1);
    forecasts = yt_window(end) * ones(12, 1);  
    forecast_matrix_3(:, i-start_forecast+1) = forecasts;
end

%% Q6
forecast_matrix_4 = zeros(12, length(quarters) - start_forecast);

for i = start_forecast:length(quarters)
    yt_window = ytrue(1:i-1);
    forecasts = zeros(12, 1);
    for h = 1:12
        forecasts(h) = mu(i) + phi^(h) * (yt_window(end) - mu(i));
    end
    forecast_matrix_4(:, i-start_forecast+1) = forecasts;
end


%% Q7
realization_matrix = NaN(12, length(quarters) - start_forecast+1);

for i = start_forecast:length(quarters)
    realization_horizon = min(length(quarters) - i + 1, 12);
    realization_matrix(1:realization_horizon, i-start_forecast+1) = ytrue(i:i+realization_horizon-1);
end

%% Q8

rmse_matrix = zeros(12, 4);  
mae_matrix = zeros(12, 4);

for i = 1:4

    if i == 1
        forecasts = forecast_matrix_1;
    elseif i == 2
        forecasts = forecast_matrix_2;
    elseif i == 3
        forecasts = forecast_matrix_3;
    else
        forecasts = forecast_matrix_4;
    end

    for h = 1:12

        true_y = realization_matrix(h, :);
        common_period = ~isnan(true_y) & ~isnan(forecasts(h, :));
        rmse_matrix(h, i) = sqrt(nanmean((forecasts(h, common_period) - true_y(common_period)).^2));
        mae_matrix(h, i) = nanmean(abs(forecasts(h, common_period) - true_y(common_period)));
    end
end


figure;
subplot(2, 1, 1);
plot(1:12, rmse_matrix(:, 1:3) ./ rmse_matrix(:, 4), 'LineWidth', 1.5);
%title('Root Mean Square Error (RMSE) - Ratio to Researcher 4');
legend('Researcher 1 - expanding window', 'Researcher 2 - rolling window', 'Researcher 3 - random walk');
ylabel('RMSE Ratio');
xlabel('Forecasting Horizon');

subplot(2, 1, 2);
plot(1:12, mae_matrix(:, 1:3) ./ mae_matrix(:, 4), 'LineWidth', 1.5);
%title('Mean Absolute Error (MAE) - Ratio to Researcher 4');
legend('Researcher 1 - expanding window', 'Researcher 2 - rolling window', 'Researcher 3 - random walk');
ylabel('MAE Ratio');
xlabel('Forecasting Horizon');

saveas(gcf, '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Plot/ps1_fig5', 'epsc');

%% Annex - Simulate 10,000 processes and compute the mean RMSE/MAE across the 4 forecast strategies.

num_iterations = 10000;
rng(8); 
all_rmse = zeros(12, 4, num_iterations);
all_mae = zeros(12, 4, num_iterations);

for iter = 1:num_iterations

    phi = 0.9;
     
    unconditional_std = 3;
    
    eps_sigma2 = (unconditional_std^2) * (1 - phi^2);
    
    
    start_year = 1948;
    end_year = 2017.75;
    quarters = (start_year:0.25:end_year)';
    
    mu = zeros(size(quarters));
    mu(quarters >= 1948 & quarters <= 1972) = 5 * (1 - phi);
    mu(quarters > 1972 & quarters <= 1995) = 0 * (1 - phi);
    mu(quarters > 1995 & quarters <= 2004) = 4.5 * (1 - phi);
    mu(quarters > 2004) = 0.5 * (1 - phi);
    
    y = zeros(size(quarters));
    epsilon = sqrt(eps_sigma2) * randn(size(quarters));
    y(1) = mu(1) + epsilon(1);
    
    for t = 2:length(quarters)
        y(t) = mu(t) + phi * y(t - 1) + epsilon(t);
    end
    
    % Q2
    
    mean_y_1948_1972 = mean(y(quarters >= 1948 & quarters <= 1972));
    mean_y_1972_1995 = mean(y(quarters > 1972 & quarters <= 1995));
    mean_y_1995_2004 = mean(y(quarters > 1995 & quarters <= 2004));
    mean_y_2004_2017 = mean(y(quarters > 2004));
        
    ytrue = y;
    
    % Q3
    start_forecast = find(quarters == 1990.00);
    
    forecast_matrix_1 = zeros(12, length(quarters) - start_forecast+1);
    
    for i = start_forecast:length(quarters)
        y_window = ytrue(1:i-2);
    
        X = [ones(length(y_window), 1), y_window];
        y = ytrue(2:i-1); 
        beta_hat = X\y;
    
        forecast_horizon = 12;
        forecasts = zeros(forecast_horizon, 1);
        for h = 1:forecast_horizon
            forecasts(h) = beta_hat(1) + beta_hat(2)^(h) * (y(end) - beta_hat(1)) ;
        end
    
         forecast_matrix_1(:, i-start_forecast+1) = forecasts;
    end
    
    % Q4
    
    forecast_matrix_2 = zeros(12, length(quarters) - start_forecast+1);    
    rolling_window_size = 60;
    
    for i = start_forecast:length(quarters)

        yt_window = ytrue(i - rolling_window_size - 2:i-2);
    
        
        X = [ones(length(yt_window), 1), yt_window];
        y = ytrue(i - rolling_window_size-1:i-1); 
        beta_hat = X\y;
    
        
        forecast_horizon = 12;
        forecasts = zeros(forecast_horizon, 1);
        for h = 1:forecast_horizon
            forecasts(h) = beta_hat(1) + beta_hat(2)^(h) * (y(end) - beta_hat(1));
        end
    
        forecast_matrix_2(:, i-start_forecast+1) = forecasts;
    end
    
    % Q5
    forecast_matrix_3 = zeros(12, length(quarters) - start_forecast+1);
    
    
    for i = start_forecast:length(quarters)
        yt_window = ytrue(1:i-1);
    
        forecasts = yt_window(end) * ones(12, 1);  
    
        forecast_matrix_3(:, i-start_forecast+1) = forecasts;
    end
    
    % Q6
    
    forecast_matrix_4 = zeros(12, length(quarters) - start_forecast);
    
    for i = start_forecast:length(quarters)
    
        yt_window = ytrue(1:i-1);
    
        forecasts = zeros(12, 1);
        for h = 1:12
            forecasts(h) = mu(i) + phi^(h) * (yt_window(end) - mu(i));
        end
    
        forecast_matrix_4(:, i-start_forecast+1) = forecasts;
    end
    
    
    % Q7
    realization_matrix = NaN(12, length(quarters) - start_forecast+1);
    
    for i = start_forecast:length(quarters)
        realization_horizon = min(length(quarters) - i + 1, 12);
        realization_matrix(1:realization_horizon, i-start_forecast+1) = ytrue(i:i+realization_horizon-1);
    end

  % Q8
    for i = 1:4
        if i == 1
            forecasts = forecast_matrix_1;
        elseif i == 2
            forecasts = forecast_matrix_2;
        elseif i == 3
            forecasts = forecast_matrix_3;
        else
            forecasts = forecast_matrix_4;
        end

        for h = 1:12
            true_y = realization_matrix(h, :);
            common_period = ~isnan(true_y) & ~isnan(forecasts(h, :));

            all_rmse(h, i, iter) = sqrt(nanmean((forecasts(h, common_period) - true_y(common_period)).^2));
            all_mae(h, i, iter) = nanmean(abs(forecasts(h, common_period) - true_y(common_period)));
        end
    end
end

mean_rmse = mean(all_rmse, 3);
mean_mae = mean(all_mae, 3);

subplot(2, 1, 1);
plot(1:12, mean_rmse(:, 1:3) ./ mean_rmse(:, 4), 'LineWidth', 1.5);
%title('Root Mean Square Error (RMSE) - Ratio to Researcher 4');
legend('Researcher 1 - expanding window', 'Researcher 2 - rolling window', 'Researcher 3 - random walk');
ylabel('RMSE Ratio');
xlabel('Forecasting Horizon');

subplot(2, 1, 2);
plot(1:12, mean_mae(:, 1:3) ./ mean_mae(:, 4), 'LineWidth', 1.5);
%title('Mean Absolute Error (MAE) - Ratio to Researcher 4');
legend('Researcher 1 - expanding window', 'Researcher 2 - rolling window', 'Researcher 3 - random walk');
ylabel('MAE Ratio');
xlabel('Forecasting Horizon');
saveas(gcf, '/Users/adrienroussetplanat/Library/CloudStorage/Dropbox/PhD/Courses/S2/Metrics_II/PS/Output/Plot/ps1_fig6', 'epsc');
