function [C_XX_shr, C_XY_shr, C_YY_shr] = shrinkPSD(R, n)

% This function shrinks the full covariance matrix until it becomes PSD,
% and returns building blocks of the resulting shrunken full covariance 
% matrix: C_XX, C_XY, C_YY.

% Anna Cichonska
% anna.cichonska@helsinki.fi




% Eigenvalues of the full covariance matrix
lambdas = eig(R); 

 
while ( min(lambdas) < 0 )    

    R = 0.999*R;                        % shrinkage
    R(logical(eye(size(R)))) = 1;       % setting diagonal elements to 1
    
    lambdas = eig(R);

end 


C_XX_shr = R(1:n, 1:n);           
C_YY_shr = R((n+1):end, (n+1):end);      
C_XY_shr = R(1:n, (n+1):end); 