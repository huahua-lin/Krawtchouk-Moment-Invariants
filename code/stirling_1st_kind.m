% This script computes Stirling numbers of the 1st kind.

function [s] = stirling_1st_kind(n, i)
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
            s = -(n-1)*stirling_1st_kind(n-1, i) + stirling_1st_kind(n-1, i-1);        
        end
    end
