function [C, area] = compute_capacitance(conn, x, y, d_elm, e0, er, u, V1, V2)
        M = size(conn,1);
        eps_e = ones(M,1)*e0; eps_e(d_elm) = er*e0;
        area  = zeros(M,1);
        Ex = zeros(M,1); Ey = zeros(M,1);
        for e = 1:M
            [~, area(e)] = assemble_S(e, x, y, conn, eps_e(e), eps_e(e), 0);
            [Ex(e), Ey(e)] = Efield(e, x, y, conn, u);
        end
        E_teo = abs(V1-V2) / (max(y)-min(y)) * 10;
        good  = d_elm(sqrt(Ex(d_elm).^2+Ey(d_elm).^2) <= 5*abs(V1-V2)/0.01);
        C = (1/(V1-V2)^2) * sum(eps_e(good).*(Ex(good).^2+Ey(good).^2).*area(good));
end