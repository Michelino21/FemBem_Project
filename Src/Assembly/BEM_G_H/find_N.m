% Calculates the normal versor N(t) = (nx(t), ny(t)) needed to evaluate the
% kernel.
%
%      met      Geometric element type (1 for linear elements, 2 for
%               quadratic).
%
%    (x1, y1) (x2, y2) (x3, y2)
%               The boundary points that identify the element. In case of
%               linear elements point (x2, y2) is the middle point between
%               (x1, y1) and (x3, y3).
%
%       t       The point in [-1, 1] where the function has to be
%               evaluated.

function [nx, ny] = find_N(met, x1, y1, x2, y2, x3, y3, t)

    % Infinitesimal increment for the numerical
    % evaluation of the normal direction.
    dt = 1e-5;

    if met == 1
        
        % Linear elements
        
        xt_p = shape_fun(met, 1, t+dt) * x1 + ...
               shape_fun(met, 3, t+dt) * x3;
        
        yt_p = shape_fun(met, 1, t+dt) * y1 + ...
               shape_fun(met, 3, t+dt) * y3;

        xt_m = shape_fun(met, 1, t-dt) * x1 + ...
               shape_fun(met, 3, t-dt) * x3;
        
        yt_m = shape_fun(met, 1, t-dt) * y1 + ...
               shape_fun(met, 3, t-dt) * y3;

    elseif met == 2
        
        % Quadratic elements
        
        xt_p = shape_fun(met, 1, t+dt) * x1 + ...
               shape_fun(met, 2, t+dt) * x2 + ...
               shape_fun(met, 3, t+dt) * x3;
           
        yt_p = shape_fun(met, 1, t+dt) * y1 + ...
               shape_fun(met, 2, t+dt) * y2 + ...
               shape_fun(met, 3, t+dt) * y3;

        xt_m = shape_fun(met, 1, t-dt) * x1 + ...
               shape_fun(met, 2, t-dt) * x2 + ...
               shape_fun(met, 3, t-dt) * x3;
           
        yt_m = shape_fun(met, 1, t-dt) * y1 + ...
               shape_fun(met, 2, t-dt) * y2 + ...
               shape_fun(met, 3, t-dt) * y3;

    else
        
        disp('Shape function error: supported element types are 1 or 2.');

    end

    ls = sqrt( (xt_p - xt_m).^2 + (yt_p - yt_m).^2 );
    
    nx = ( yt_p - yt_m ) ./ ls;
    ny = ( xt_m - xt_p ) ./ ls;

end