function S_XY = normalize_Sxy(S_XY_raw, se, N)

% This function normalizes regression coefficients 
% using their standard errors 'se' and the number of samples 'N'.

% Anna Cichonska
% anna.cichonska@helsinki.fi

% The software is for academic purposes only.
% Commercial use is not allowed.
% The software is provided "as is", without warranty of any kind.



S_XY =  (1/sqrt(N)) .* (S_XY_raw./se);

