% Calculates the radius R(t) = (rx(t), ry(t)) needed to evaluate the
% kernel.
%
%      met      Geometric element type (1 for linear elements, 2 for
%               quadratic).
%
%    (xi, yi)   The reference point.
%
%    (x1, y1) (x2, y2) (x3, y2)
%               The boundary points that identify the element. In case of
%               linear elements point (x2, y2) is the middle point between
%               (x1, y1) and (x3, y3).
%
%       t       The point in [-1, 1] where the function has to be
%               evaluated.

function [rx, ry] = find_R(met, xi, yi, x1, y1, x2, y2, x3, y3, t)

    if met==1

        % Linear elements

        xt = shape_fun(met, 1, t) * x1 + ...
             shape_fun(met, 3, t) * x3;
        
        yt = shape_fun(met, 1, t) * y1 + ...
             shape_fun(met, 3, t) * y3;

    elseif met == 2

        % Quadratic elements

        xt = shape_fun(met, 1, t) * x1 + ...
             shape_fun(met, 2, t) * x2 + ...
             shape_fun(met, 3, t) * x3;
        
        yt = shape_fun(met, 1, t) * y1 + ...
             shape_fun(met, 2, t) * y2 + ...
             shape_fun(met, 3, t) * y3;
    
    else
        
        disp('Shape function error: supported element types are 1 or 2.');

    end

    rx = xt - xi;
    ry = yt - yi;

end