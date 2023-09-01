% prepStep.m -- This file is part of 3DKMI.
%
% This script follows the processes of calculating the 3D weight function, the 3D norms, and the coefficients.

function const = prepStep(N, order)
    const.order = order;

    % Recursively compute the square-root of the 1D weight function w(x;p,N-1) 
    p = 0.5;
    w = zeros(N,1);
    w(1) = 1;
    for x = 0:floor((N-1)/2)-1
        w(x+2) = w(x+1) *  sqrt((N-1-x) / (x+1)) ;
    end
    w = w * p^((N-1)/2);

    % Use symmetry to produce the second half of w(x;p,N-1) from its first half
    w( N-floor(N/2)+1 : 1 : N ) = w( floor(N/2) : -1 : 1 );
    
    % Compute the 3D weight function for px = py = pz = 0.5
    wlen = size(w,1);
    ww = w*w';
    const.Wc = zeros(wlen, wlen, wlen);
    for i = 1:wlen
        const.Wc(:,:,i) = ww(:,:) * w(i);
    end

    % Compute the norms rho(n;p,N-1)
    rho = zeros(order+1,1);
    for n = 0:order
        rho(n+1) = (-1)^n * ((1-p) / p)^n * factorial(n) / pochhammer(-(N-1),n);
    end
    rholen = size(rho,1);
    rhorho = rho * rho';
    const.rho = zeros(rholen, rholen, rholen);
    for i = 1:rholen
        const.rho(:,:,i) = rhorho(:,:) * rho(i);
    end

    % Compute the coefficients a(i,n,p,N-1)
    const.a = zeros(order+1,order+1);
    for n = 0:order
        for i = 0:n
            for k = i:n
%                 const.a(i+1,n+1) = const.a(i+1,n+1) + pochhammer(-n,k)/(pochhammer(-N+1,k) * factorial(k)) * (-2)^k * stirling_1st_kind(k,i);
                const.a(i+1,n+1) = const.a(i+1,n+1) + ((-1)^k * factorial(n) * factorial(N-1-k))/(p^k * factorial(N-1) * factorial(n-k) * factorial(k)) * stirling_1st_kind(k,i);
            end
        end
    end