% It has the same goal of findG(), but it handles the logarithmic
% singularity in case of linear elements, which can be integrated
% analitically.

function [g] = findG_lin_s(x1, y1, x2, y2)

    lm = sqrt( (x2 - x1)^2 + (y2 - y1)^2 );

    g = lm/(4*pi) * (log(lm) - 1.5);

end
