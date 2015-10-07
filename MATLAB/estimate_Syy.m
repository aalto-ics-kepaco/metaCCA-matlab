function S_YY = estimate_Syy(S_XY, N, ind)

% This function 
% 1) normalizes regression coefficients ('normalize_Sxy') if the univariate 
%    analsysis has been performed on non-standardised data;
% 2) uses Pearson correlation to compute phenotypic correlation matrix S_YY
%    based on univariate summary statistics S_XY.

% Anna Cichonska
% anna.cichonska@helsinki.fi



% Validating if trait ids are provided as strings
if size(S_XY.textdata,2)-3 ~= size(S_XY.data,2)
    error('Trait ids in the header line should be given as strings.');
end

trait_ids_betas = S_XY.textdata(1,4:2:end);
trait_ids_se    = S_XY.textdata(1,5:2:end);

% Validating if trait ids of the correspondung regression coefficients 
% and standard errors match
if isequal(trait_ids_betas, trait_ids_se) == 0
    error('Trait ids of regression coefficients and standard errors do not match');
end


% Validating if trait ids are unique
if length(unique(trait_ids_betas)) ~= length(trait_ids_betas)
   error('Trait ids are not unique!');
end



Betas = S_XY.data(:, 1:2:end);
SE    = S_XY.data(:, 2:2:end);


if ind == 0
    S_XY_norm = normalize_Sxy(Betas, SE, N);
elseif ind == 1
    S_XY_norm = Betas;
else
    error('Wrong indicator of the data type. Please select >>0<< or >>1<<. ')
end


S_YY = corr(S_XY_norm);

S_YY = [trait_ids_betas', num2cell(S_YY)];
