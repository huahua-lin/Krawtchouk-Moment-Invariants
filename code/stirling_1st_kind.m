function [ s ] = stirling_1st_kind( n, i )
% This function computes Stirling numbers of the 1st kind, as described by
% Wikipedia.
%
% Input:
%   n = argument 1
%   i = argument 2
%
% Output:
%   s = s(n,i), where s is the Stirling number of the first kind, defined
%       by the recurrence relation:
%           s(n, i) = -(n-1)*s(n-1 , i) + s(n-1, i-1)
%       with the initial conditions
%           s(0,0) = 1,   s(n,0) = s(0,i) = 0, if n > 0
%       See the Wikipedia article: https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind
if n == 0
    if i == 0
        s = 1;
    else %i > 0
        s = 0;
    end
else %n > 0
    if i == 0
       s = 0; 
    else %i > 0
        %this is the general case ==> use the recursive formula
        s = -(n-1)*stirling_1st_kind( n-1, i ) + stirling_1st_kind( n-1, i-1 );        
    end
end
