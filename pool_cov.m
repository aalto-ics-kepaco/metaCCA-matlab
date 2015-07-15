function [C_XX, C_YY, C_XY, sum_N] = pool_cov(varargin)

% This function pools covariance matrices of the same type;
% 'varargin' argument is a cell array that contains the inputs, where each input is in its own cell
% There are 4 types of inputs, in order: C_XX, C_YY, C_XY, N

% Anna Cichonska
% anna.cichonska@helsinki.fi


nr_in = nargin;          % number of inputs
nr_s  = nargin/4;         % number of studies


counter_XX = 0;
counter_YY = 0;
counter_XY = 0;
denom      = 0;
sum_N      = 0;


for i = 1:nr_s
    counter_XX = counter_XX + ((varargin{nr_in-nr_s+i}-1) * varargin{i});               
    counter_YY = counter_YY + ((varargin{nr_in-nr_s+i}-1) * varargin{i+nr_s});
    counter_XY = counter_XY + ((varargin{nr_in-nr_s+i}-1) * varargin{i+2*nr_s});                              
    
    denom = denom + (varargin{nr_in-nr_s+i}-1);                                         
    
    sum_N = sum_N + varargin{nr_in-nr_s+i};
end


C_XX = counter_XX/denom;
C_YY = counter_YY/denom;
C_XY = counter_XY/denom;
