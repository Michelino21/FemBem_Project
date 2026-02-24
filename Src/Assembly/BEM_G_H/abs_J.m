% Calculates the jacobian for a given element (the element type is
% automatically recognized by the given point-coordinates).
%
%        t      The point in [-1, 1] where the function has to be
%               evaluated.
%
%    (x1, y1) (x2, y2) (x3, y2)
%               The boundary points that identify the element. In case of
%               linear elements point (x2, y2) is the middle point between
%               (x1, y1) and (x3, y3).

function [J] = abs_J(x1, y1, x2, y2, x3, y3, t)

    b1 = (x3 - x1) / 2;
    b2 = (x1 - 2*x2 + x3) / 2;
    c1 = (y3 - y1) / 2;
    c2 = (y1 - 2*y2 + y3) / 2;

    J = sqrt((b1 + 2*b2*t).^2 + (c1 + 2*c2*t).^2);

end