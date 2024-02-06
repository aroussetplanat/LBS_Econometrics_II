%% 1.Replicating Kilian (2009)
 %1.1 
dt = load('../Data/Kilian_Data_Updated.mat');
dt.dates= datetime(dt.dates, 'ConvertFrom', 'datenum');
Shocknames={'Oil Supply Shock','Aggregate Demand Shock','Oil-Specific Demand Shock'};
[T,n] = size(dt.data);

f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])


for i = 1:n
    subplot(numel(dt.varNames), 1, i);
    plot(dt.dates, dt.data(:, i),'b-');   
    title(dt.varNames{i});

end
saveas(gcf, '../Output/Plot/ps2_fig1', 'epsc');
%% 1.2
constant=1;
p=12;
results=VAR_OLS(dt.data,constant,p); %dt.data
%% 1.3
A0 = chol(results.Sigma, 'lower');

%% 1.4
hmax=T-p;
hmaxtoplot=24;
horizon=0:hmax;

cumulateWhich=1;
IRF=getIRF(results.B, A0,constant, hmax, cumulateWhich);

f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

x  = 1;

for jj = 1:n
 for ii = 1:n
    subplot(3,3,x)
    
    if ismember(ii, [1, 2, 3])  && jj==1 
    y = -squeeze(IRF(ii,jj,:));
    else
    y = squeeze(IRF(ii,jj,:));
    end 

    plot(horizon, y, 'b-', horizon, zeros(size(horizon)),'k--')
    
    xlim([0,hmaxtoplot])
    if ii == 1
    ylim([-2.5,1.5]) 
    ylabel('Percent')
    else
    ylim([-5,15])
    end
    
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 8); 
    set(gca,'XTick',[0:round(hmaxtoplot/4.33):hmaxtoplot]')
    set(gca,'XTickLabel',num2str([0:round(hmaxtoplot/4.33):hmaxtoplot]'))    
    title(strcat(dt.varNames{ii},{' '},'to',{' '},Shocknames{jj},' Shock'))
    set(gca,'Layer','top')

    if x > 6
            xlabel('Months')
    end
    x  = x+1;
    
 end
    
end

saveas(gcf, '../Output/Plot/ps2_fig2', 'epsc');

%% 1.5
structural_u=getShocks(dt.data,results.B,A0,constant, p);

f = figure;
clf
figSize = [12 6];
for i = 1:n
    subplot(n, 1, i);
    plot(dt.dates(p+1:end), structural_u(:, i),'b-',dt.dates(p+1:end),zeros(size(structural_u(:, i))), 'k--');   
    title(Shocknames{i});
end
saveas(gcf, '../Output/Plot/ps2_fig3', 'epsc');

u_yearly = retime(timetable(dt.dates(p+1:end), structural_u), 'yearly', 'mean');

f = figure;
clf
figSize = [12 6];
for i = 1:n
    subplot(n, 1, i);
    plot(u_yearly.Time, u_yearly{:, 1}(:,i),'b-',u_yearly.Time, zeros(size(u_yearly.Time)), 'k--');   
    title(Shocknames{i});
end
      
saveas(gcf, '../Output/Plot/ps2_fig4', 'epsc');


%% 1.6
FEVD=getFEVD(IRF,hmaxtoplot);
for i = 1:n
        figure;
        clf
        bar(0:hmaxtoplot, squeeze(FEVD(i, :, :)),'stacked');
        xlabel('Forecast Horizon');
        ylabel('FEVD');
        legend(dt.varNames, 'Location', 'best'); 
        %title([dt.varNames{i}]);
        axis([0 25 0 1]);
        xlim([-0.5, hmaxtoplot+0.5]);
        saveas(gcf, ['../Output/Plot/ps2_fig', num2str(i+4)], 'epsc');

end


%% 1.7
whichVar=3;
HD_store=nan(size(structural_u));
dt.data_demean=dt.data-mean(dt.data);

for t = 1:T-p
    HD_store(t,:) = getHDs(IRF(:,:,1:t),n,structural_u(1:t,:),whichVar);
end 

sum_HD_store=sum(HD_store,2);

f = figure;
clf
figSize = [12 6];
set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

subplot(2, 1, 1);
plot(dt.dates(p+1:end), HD_store(:, 1), 'b-', ...
     dt.dates(p+1:end), HD_store(:, 2), 'r-', ...
     dt.dates(p+1:end), HD_store(:, 3), 'g-', ...
     dt.dates(p+1:end), dt.data_demean(p+1:end, whichVar), 'k-', ...
     'LineWidth', 2); 
hold on; 
line([dt.dates(p+1), dt.dates(end)], [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); 
legend(Shocknames{1}, Shocknames{2}, Shocknames{3},  [dt.varNames{whichVar}, ' - net of deterministic component']);
subplot(2, 1, 2);

plot(dt.dates(p+1:end), dt.data_demean(p+1:end, whichVar), 'k-', ...
     dt.dates(p+1:end), sum_HD_store(:, 1), 'c-', ...
     dt.dates(p+1:end), dt.data_demean(p+1:end, whichVar), 'k-', ...
     'LineWidth', 2);
hold on; 
line([dt.dates(p+1), dt.dates(end)], [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); 
     legend('Real Oil Price - net of deterministic component'  , 'Sum of the contribution of the three structural shocks'); 

saveas(gcf, '../Output/Plot/ps2_fig8', 'epsc');


%% 1.8
M=2000;
[T, n] = size(dt.data);
IRF_bands = zeros(n, n, hmaxtoplot+1,M);
Yr = zeros(T, n);
Yr(1:p, :) = dt.data(1:p, :);
if constant == 1
    c = results.B(1, :)';    
elseif constant == 0
    c=zeros(n,1);
end

for m = 1:M
    Yr(p+1:end, :) = 0;
    idx = randsample(T-p, T-p, true);
    V = results.residuals(idx, :);
    for j=p+1:T   
        for k = 1:p                      
            Yr(j,:) = Yr(j,:) + Yr(j-k, :)*results.C{k}';
        end
        Yr(j,:) = c' + Yr(j,:) + V(j-p, :);
    end 
    results_m=VAR_OLS(Yr,constant,p);
    A0_m = chol(results_m.Sigma, 'lower');
    IRF_m=getIRF(results_m.B, A0_m,constant, hmaxtoplot, cumulateWhich);
    IRF_bands(:,:,:,m) = IRF_m;

end 

bands  = [5,50,95];
f = figure;
clf
figSize = [12 6];

set(f, 'PaperUnits', 'inches');
set(f, 'Units','inches');
set(f, 'PaperSize', figSize);
set(f, 'PaperPositionMode', 'auto');
set(f, 'Position', [0 0 figSize(1) figSize(2)])

x  = 1;

for jj = 1:n
 for ii = 1:n
    subplot(3,3,x)

    if ismember(ii, [1, 2, 3])  && jj==1 
    percentiles = prctile(-squeeze(IRF_bands(ii,jj,1:hmaxtoplot+1,:)),bands,2);
    else
    percentiles = prctile(squeeze(IRF_bands(ii,jj,1:hmaxtoplot+1,:)),bands,2);
    end 
        
    plotConfidenceBandsBlue(0:hmaxtoplot,percentiles,'b');
    hold on
    plot(0:hmaxtoplot,zeros(hmaxtoplot+1,1),'--k')

    xlim([0,hmaxtoplot])
    if ii == 1
    ylim([-2.5,1.5])  
    ylabel('Percent')
    else
    ylim([-5,15])
    end

    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 8); 
    set(gca,'XTick',[0:round(hmaxtoplot/4.33):hmaxtoplot]')
    set(gca,'XTickLabel',num2str([0:round(hmaxtoplot/4.33):hmaxtoplot]'))    
    box on
    title(strcat(dt.varNames{ii},{' '},'to',{' '},Shocknames{jj},' Shock'))
    set(gca,'Layer','top')

    if x > 6
            xlabel('Months')
    end
    x  = x+1;
    
 end
    
end
saveas(gcf, '../Output/Plot/ps2_fig9', 'epsc');

