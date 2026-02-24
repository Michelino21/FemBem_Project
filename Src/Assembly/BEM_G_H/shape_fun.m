% Calculates the value of the shape function in the given point t, which
% has to be in [-1, 1].
%
%    met   Geometric element type (1 for linear elements, 2 for quadratic).
%     c    Shape function type (2 is allowed only with quadratic elements).
%     t    The point in [-1, 1] where the function has to be evaluated.

function [sf] = shape_fun(met, c, t)

    if met == 1
        
        % Linear elements
        
        if c == 1
            sf = 0.5 * (1-t);
        elseif c == 3
            sf = 0.5 * (1+t);
        else
            disp('Shape function error: c cannot be 2')
        end
        
    elseif met == 2
        
        % Quadratic elements
        
        if c == 1
            sf = -0.5 * t .* (1-t);
        elseif c == 2
            sf = (1-t) .* (1+t);
        elseif c==3
            sf = +0.5 * t .* (1+t);
        end
     
    else
        
        disp('Shape function error: supported element types are 1 or 2.');
        
    end

end
    
        
            