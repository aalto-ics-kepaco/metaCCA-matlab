function [r, a, b, wilks, chi, df, p_val] = my_cca(C_XX, C_YY, C_XY, N) 

% This function performs Canonical Correlation Analysis (CCA) based on the
% building blocks of the full covariance matrix (C_XX, C_YY, C_XY), as well as
% significance testing of the canonical correlations.

% Anna Cichonska
% anna.cichonska@helsinki.fi



K = C_XX^(-0.5) * C_XY * C_YY^(-0.5);

% SVD Singular Value Decomposition
[U,S,V] = svd(K);        % K = U*S*V'
% diag(S) - r (canonical correlation)


min_r = min(size(C_XX,1), size(C_YY,1));    
a = nan(size(C_XX,1), min_r);   % canonical weights a and b
b = nan(size(C_YY,1), min_r);
r = nan(1, min_r);              % canonical correlation
wilks = nan(1,min_r);           % Wilk's lambda
lambdas = S.^2;                 % r^2 = eigenvalue corresponding to r
df = nan(1,min_r);              % degrees of freedom

p = size(C_XX,1);
q = size(C_YY,1);
k=0;

for i = 1:min_r
    a(:,i) = C_XX^(-0.5) * U(:,i);
    b(:,i) = C_YY^(-0.5) * V(:,i);
    r(1,i) = (a(:,i)'*C_XY*b(:,i)) / (sqrt(a(:,i)'*C_XX*a(:,i)) * sqrt(b(:,i)'*C_YY*b(:,i)));    
    

    % Wilk's lambda
    prod = 1;
    for j = i:min_r
        prod_temp = 1-lambdas(j,j);
        prod = prod*prod_temp;
    end
    wilks(i) = prod;
    df(i) = (p-k)*(q-k);       % df for Bartlett Chi-square
    k = k+1;
end


    
% Bartlett Chi-square approximation to Wilk's Lambda
chi = -((N-1)- 0.5*(p+q+1)) * log(wilks);

% H0: all canonical correlations = 0
% p-value
 p_val = chi2cdf_my(chi,df);

