function [ Ylag ] = mlag2 ( Y_raw , p )
% Mlag2 creates a matrix of p lags from Y_raw
Y = Y_raw ;
[T , N ]= size ( Y );
Ylag = zeros (T , N * p );
for i =1: p
Ylag ( p +1: T ,( N *( i -1)+1): N * i )= Y ( p +1 - i :T -i ,1: N );
end