function opt_threshold = optimal_threshold(d_full,R,n,i)

% This function finds an appropriate shrinkage level of the full 
% covariance matrix using an elbow criterion.

% Anna Cichonska
% anna.cichonska@helsinki.fi




% Monitoring the behaviour of the leading canonical correlation
% value during shrinkage of the full covariance matrix
while  abs(diff(d_full(end-1:end)))/d_full(end-1) > 0.005   
    
    i = i+1;
    R = 0.999*R;                        
    R(logical(eye(size(R)))) = 1;       
    
    C_XX_new_temp = R(1:n, 1:n);           
    C_YY_new_temp = R((n+1):end, (n+1):end);      
    C_XY_new_temp = R(1:n, (n+1):end); 
                     
    K = C_XX_new_temp^(-1/2) * C_XY_new_temp * C_YY_new_temp^(-1/2);
    d_full(i+1) = max(svd(K)); 
end


% Finding an appropriate amount of shrinkage
d_temp = abs(diff(d_full)); d_temp = d_temp(3:end);  d_temp = d_temp./d_full(3:(end-1));
if length(d_temp)<2  
    opt_threshold = 0.005;
else
    x = 2:(length(d_temp)+1);
    y = d_temp;

    % Normalizing x and y to [0,1]
    x_n = (x-min(x))/(max(x)-min(x));
    y_n = (y-min(y))/(max(y)-min(y));

    coordinates = [x_n;y_n]';

    % figure; plot(x_n,y_n,'o')

    % Computing Euclidean distance of each point to the origin of the plot
    % of the percent change of canonical correlations between subsequent shrinkage 
    % iterations versus iteration number
    for coord_i = 1:length(coordinates)
        dist_to_origin(coord_i) = norm( coordinates(coord_i,:) - [0,0]);
    end

    % Taking the point which has the smallest distance to the origin
    % of the plot
    [~, id_min_dist] = min(dist_to_origin);

    opt_threshold = y(id_min_dist);
end