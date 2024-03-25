function out   = TransformSeries(rawdata,tcode)
% =========================================================================
% Transforms raw data based on each series' transformation code.
% =========================================================================

[T,N] = size(rawdata);
out = NaN(T,N);
% Value close to zero to 
small = 1e-6;

for i = 1:N % loop over each variables 
    x=rawdata(:, i);
    % Number of observations (including missing values)
    n = size(x, 1);
    y=NaN*ones(n,1);

    % Perform transformation based on transformation code
    if tcode(i) == 1 % Level (no transformation)
        y = x;
    
    elseif tcode(i) == 2 % First difference
        y(2:n)=x(2:n,1)-x(1:n-1,1);
    
    elseif tcode(i) == 3 % Second difference
        y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
    
    elseif tcode(i) == 4 % Natural log
        if min(x) < small
            y=NaN; 
        else
            y=log(x);
        end
    
    elseif tcode(i) == 5 % First difference of natural log
        if min(x) > small
            x=log(x);
            y(2:n)=x(2:n)-x(1:n-1);
        end
    
    elseif tcode(i) == 6 % Second difference of natural log
        if min(x) > small
            x=log(x);
            y(3:n)=x(3:n)-2*x(2:n-1)+x(1:n-2);
        end
    
    elseif tcode(i) == 7 % First difference of percent change
        y1(2:n)=(x(2:n)-x(1:n-1))./x(1:n-1);
        y(3:n)=y1(3:n)-y1(2:n-1);
    end
    
        % Append transformed series to output
        out(:,i)=y;
end

end

