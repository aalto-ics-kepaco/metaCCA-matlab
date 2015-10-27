function S_XY = normalize_Sxy(S_XY_raw, se, N)

% This function normalizes regression coefficients 
% using their standard errors 'se' and the number of samples 'N'.

% Anna Cichonska
% anna.cichonska@helsinki.fi



S_XY =  (1/sqrt(N)) .* (S_XY_raw./se);

