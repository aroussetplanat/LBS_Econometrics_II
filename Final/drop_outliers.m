function Y=drop_outliers(X)
% =========================================================================
% replaces outliers with the value NaN as in McCracken and Ng (2020).
% A data point x of a variable is considered an outlier 
% if abs(x-median)>10*interquartile_range.
% =========================================================================

% Compute median and quartiles of each variable
median_X = nanmedian(X);
Q = prctile(X, [25, 75]);

% Calculate interquartile range (IQR) of each series
IQR = Q(2, :) - Q(1, :);

% Determine outliers 
outlier = abs(X - median_X) > (10 * IQR);

% Replace outliers with NaN
Y = X;
Y(outlier) = NaN;
end

