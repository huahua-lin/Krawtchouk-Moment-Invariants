function invariants = extractFeatures(fs, xs, ys, zs, const)    
    N = size(const.Wc,1);
    
    % Translate const.Wc computed in prepStep so its geometric center moves from 
    % the unique point ((N-1)/2,(N-1)/2,(N-1)/2) to the point-of-interest (xs,ys,zs).
    W = zeros(N,N,N);
    xmove = round( xs - (N-1)/2 );
    if xmove > 0
        W( 1+xmove : N, : , : ) = const.Wc( 1 : N-xmove, : , : ); % 1~S -> (1+xmove)~(S+xmove), outside the grid S -> 0
    elseif xmove < 0
        W( 1 : N+xmove, : , : ) = const.Wc( 1-xmove : N, : , : );
    else
        W = const.Wc;
    end
    
    Wt = zeros(N,N,N);
    ymove = round( ys - (N-1)/2 );
    if ymove > 0
        Wt( : , 1+ymove : N, : ) = W( : , 1 : N-ymove, : );
    elseif ymove < 0
        Wt( : , 1 : N+ymove, : ) = W( : , 1-ymove : N, : );
    else
        Wt = W;
    end
    clear W;

    Ws = zeros(N,N,N);
    zmove = round(zs - (N-1)/2);
    if zmove > 0
        Ws(:,:,1+zmove:N) = Wt(:,:,1:N-zmove);
    elseif zmove < 0
        Ws(:,:,1 : N+zmove) = Wt(:,:,1-zmove : N);
    else
        Ws=Wt;
    end
    clear Wt;

    % weighted image ftilde
    ftilde = fs .* const.Wc;
    clear fs;

    x = 0:1:N-1;
    y = 0:1:N-1;
    z = 0:1:N-1;
    sumftilde12 = squeeze(sum(sum(ftilde,1), 2));
    sumftilde23 = sum(sum(ftilde,2), 3);
    sumftilde13 = sum(sum(ftilde,3), 1);

    % geometric moments of ftilde (5)
    M000 = sum(sumftilde23);
    M100 = sum(x'.*sumftilde23);
    M010 = sum(y .*sumftilde13);
    M001 = sum(z' .*sumftilde12);

    % the center of mass (xtilde, ytilde, ztilde) of ftilde
    xtilde = M100 / M000;
    ytilde = M010 / M000;
    ztilde = M001 / M000;
    
    xMinuSxtilde = x - xtilde;
    clear x;
    yMinuSytilde = y - ytilde;
    clear y;
    zMinuSztilde = z - ztilde;
    clear z;
    
    sumftilde1 = squeeze(sum(ftilde, 1));
    sumftilde2 = squeeze(sum(ftilde, 2));
    sumftilde3 = sum(ftilde, 3); 

    % central moments (6)
    mu200 = sum((xMinuSxtilde.^2)' .* sumftilde23);
    mu020 = sum((yMinuSytilde.^2) .* sumftilde13);
    mu002 = sum((zMinuSztilde.^2)' .* sumftilde12);
    mu110 = sum(sum(xMinuSxtilde' * yMinuSytilde .* sumftilde3));
    mu101 = sum(sum(xMinuSxtilde' * zMinuSztilde .* sumftilde2));
    mu011 = sum(sum(yMinuSytilde' * zMinuSztilde .* sumftilde1));

    % construct the inertia matrix (8)
    inertiaXX = mu020 + mu002;
    inertiaYY = mu200 + mu002;
    inertiaZZ = mu200 + mu020;
    inertiaXY = -mu110;
    inertiaXZ = -mu101;
    inertiaYZ = -mu011;
    inertia = [inertiaXX inertiaXY inertiaXZ; inertiaXY inertiaYY inertiaYZ; inertiaXZ inertiaYZ inertiaZZ];
    clear inertiaXX inertiaYY inertiaZZ inertiaXY inertiaXZ inertiaYZ;
    
    % obtain the unique eigenvectors
    [U, ~] = eig(inertia);
    centroid = [xtilde, ytilde, ztilde];
    if centroid * U(:,1) < 0
        U(:,1) = -U(:,1);
    end
    if centroid * U(:,2) < 0
        U(:,2) = -U(:,2);
    end
    if centroid * U(:,3) < 0
        U(:,3) = -U(:,3);
    end

    X = repmat(xMinuSxtilde', [1 length(yMinuSytilde) length(zMinuSztilde)]);
    Y = repmat(yMinuSytilde, [length(xMinuSxtilde) 1 length(zMinuSztilde)]);
    Z = permute(repmat(zMinuSztilde, [length(xMinuSxtilde) 1 length(zMinuSztilde)]), [3 1 2]);

    % the eigenvectors define the rows of the rotation matrix
    phi1 = U(1,1)*X + U(2,1)*Y + U(3,1)*Z;
    phi2 = U(1,2)*X + U(2,2)*Y + U(3,2)*Z;
    phi3 = U(1,3)*X + U(2,3)*Y + U(3,3)*Z;

    % geometric moment invariants (10)
    V = zeros(const.order+1,const.order+1,const.order+1);
    ftildeAi = ftilde;
    for i = 0:const.order
        if i ~= 0
            ftildeAi = phi1.*ftildeAi;
        end
        ftildeAiBj = ftildeAi;
        for j = 0:const.order
            if j ~= 0
                ftildeAiBj = phi2.*ftildeAiBj;
            end
            ftildeAiBjCk = ftildeAiBj;
            for k = 0:const.order
                if i+j+k <= const.order
                    if k ~= 0
                        ftildeAiBjCk = phi3.*ftildeAiBjCk;
                    end
                    V(i+1,j+1,k+1) = sum(sum(sum(ftildeAiBjCk)));
                end
            end
        end
    end
    clear phi1 phi2 phi3 ftilde ftildeAi ftildeAiBj ftildeAiBjCk;
    V = V ./ M000;
    
    % 3D weighted Krawtchouk moment invariants
    Q = zeros(const.order+1,const.order+1,const.order+1);
    num = length(allVL1(const.order,const.order,'<='));
    invariants = zeros(1, num - 7);
    idx = 1;
    for n = 0:const.order
        for m = 0:const.order
            for l = 0:const.order
                if n+m+l <= const.order
                    for i = 0:n
                        for j = 0:m
                            for k = 0:l
                                Q(n+1,m+1,l+1) = Q(n+1,m+1,l+1) + const.a(i+1,n+1) * const.a(j+1,m+1) * const.a(k+1,l+1) * V(i+1,j+1,k+1);
                            end
                        end
                    end
                    Q(n+1,m+1,l+1) = Q(n+1,m+1,l+1) ./ sqrt(const.rho(n+1,m+1,l+1));
                    if ~((n==0&&m==0&&l==0) || (n==1&&m==0&&l==0) || (n==0&&m==1&&l==0) || (n==0&&m==0&&l==1) || (n==0&&m==1&&l==1) || (n==1&&m==0&&l==1) || (n==1&&m==1&&l==0))
                        invariants(1,idx) = Q(n+1,m+1,l+1);
                        idx = idx + 1;
                    end
                end
            end
        end
    end
    clear V;


    