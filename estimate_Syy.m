function S_YY = estimate_Syy(S_XY)

% This function uses Pearson correlation to compute
% phenotypic correlation matrix S_YY based on summary statistics S_XY.

% Anna Cichonska
% anna.cichonska@helsinki.fi

S_YY = corr(S_XY);




 