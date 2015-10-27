function [C_XX_shr, C_XY_shr, C_YY_shr] = shrink_plus(R, n)

% This function shrinks the full covariance matrix beyond the level 
% guaranteeing its PSD property.
% It returns building blocks of the resulting shrunken full covariance 
% matrix: C_XX, C_XY, C_YY.

% Anna Cichonska
% anna.cichonska@helsinki.fi



% Eigenvalues of the full covariance matrix
lambdas = eig(R); 

% Building blocks of the full covariance matrix
C_XX = R(1:n, 1:n);           
C_YY = R((n+1):end, (n+1):end);      
C_XY = R(1:n, (n+1):end); 


d(1) = 10^10;
K = C_XX^(-1/2) * C_XY * C_YY^(-1/2);   
d(2) = max(svd(K)); 
i = 1;

% Finding an appropriate shrinkage level
opt_threshold = optimal_threshold(d,R,n,i);


% Shrinkage

while (  (abs(d(i)-d(i+1)))/d(i)  > opt_threshold   ||  (min(lambdas)<0 ) )    

    i = i+1;
    R = 0.999*R;                        % shrinkage
    R(logical(eye(size(R)))) = 1;       % setting diagonal elements to 1
    lambdas = eig(R);

    C_XX_new_gen = R(1:n, 1:n);           
    C_YY_new_gen = R((n+1):end, (n+1):end);      
    C_XY_new_gen = R(1:n, (n+1):end); 

    K = C_XX_new_gen^(-1/2) * C_XY_new_gen * C_YY_new_gen^(-1/2);
    d(i+1) = max(svd(K)); 

end 


C_XX_shr = R(1:n, 1:n);           
C_YY_shr = R((n+1):end, (n+1):end);      
C_XY_shr = R(1:n, (n+1):end); 