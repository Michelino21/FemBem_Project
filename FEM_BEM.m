%% FEM_BEM.m
% Coupling FEM-BEM per condensatore a facce piane parallele 2D
clear all; close all; clc;
addpath(genpath('src'));

SALVATAGGIO = 1;

%% SALVATAGGI
if(SALVATAGGIO == 1)

    % Scegli il nome della cartella
    folder_name = './Documentation/figures/coarse';  % cambia qui per ogni run
    
    % Crea la cartella se non esiste
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end
    
    % Nome del file di log
    log_file = fullfile(folder_name, 'console_output.txt');

    % se esiste, cancellalo
    if exist(log_file, 'file')
        delete(log_file);
    end
    
    % Attiva il salvataggio della console
    diary(log_file)
    diary on

end
%% === PROBLEM PARAMETERS ===
dx   = 0.10;          % lunghezza piastre
dy   = 0.01;          % separazione piastre  →  rapporto dx/dy = 10
delx = 0.005;         % passo mesh (20 elementi lungo dx, 2 lungo dy)
dely = delx;

margine = 3 * dy;     % = 0.03 m per lato
Lx = dx + 2*margine;  % = 0.16 m
Ly = dy + 2*margine;  % = 0.07 m

e0  = 8.85e-12;
er  = 1;
V1  = 1;   % potenziale piastra superiore
V2  = 0;   % potenziale piastra inferiore

%% === SOLVER PARAMETERS ===
elm_type = 1; % 1: triangular, 2: quadrilateral %% Only triangular

%% AVVIO TIMER
tic
%% === MESH ===
if elm_type == 1  % triangular
    [conn, x, y, xmid, ymid, outb_node, ptop_node, pbottom_node, d_elm] = ...
    meshgen_cap(Lx, Ly, dx, dy, delx, dely);    
else
   disp('Element type can only be triangular.');
   return;
end

M = size(conn, 1);
N = length(x);

eps_e = ones(M,1) * e0;
eps_e(d_elm) = er * e0;

%% === PARTIZIONE DEI NODI ===
D_nodes = [ptop_node;    pbottom_node];
B_nodes = outb_node;
F_nodes = setdiff(setdiff(1:N, D_nodes), B_nodes)';

%% === ASSEMBLAGGIO S (matrice FEM globale) ===
S = sparse(N, N);
area = zeros(M,1);

for e = 1:M
    if elm_type == 1  % triangular
        Ne = 3;
        [Ae, area(e)] = assemble_S(e, x, y, conn, eps_e(e), eps_e(e), 0);
    end
    
    for i = 1:Ne
        ig = conn(e,i);
        for j = 1:Ne
            jg = conn(e,j);
            S(ig,jg) = S(ig,jg) + Ae(i,j);
        end
    end
end

%% === BLOCCHI DI S ===
S_FF = S(F_nodes, F_nodes);
S_FB = S(F_nodes, B_nodes);
S_FD = S(F_nodes, D_nodes); %RHS
S_BF = S(B_nodes, F_nodes);
S_BB = S(B_nodes, B_nodes);
S_BD = S(B_nodes, D_nodes); %RHS

%% === RHS ===

u_D = [V1 * ones(length(ptop_node),    1);
       V2 * ones(length(pbottom_node),  1)];

%% === ASSEMBLAGGIO T_BB ===
T_BB = assemble_T(B_nodes, x, y);

%% === BEM ===
% Riordina i nodi di Gamma come percorso chiuso antiorario
B_ordered = reorder_boundary(B_nodes, x, y);

% Estrai coordinate nel giusto ordine per il BEM
xB = x(B_ordered);
yB = y(B_ordered);

[G_bem, H_bem] = assemble_G_H(1, xB, yB); %1 means linear element

% Riporta G e H nell'ordine di B_nodes (non B_ordered)
% costruisce la permutazione
[~, perm] = ismember(B_ordered, B_nodes);
G_bem = G_bem(perm, perm);
H_bem = H_bem(perm, perm);

%% === SISTEMA ACCOPPIATO ===
% Nota: segno + perché n_FEM = n_BEM
K = [S_FF,                        S_FB; ...
     S_BF,   S_BB - T_BB * (G_bem \ H_bem)];

rhs = -[S_FD; S_BD] * u_D;

%% === SOLUZIONE ===
sol   = K \ rhs;
u_F   = sol(1:length(F_nodes));
u_B   = sol(length(F_nodes)+1:end);

% Ricostruzione del potenziale globale
u = zeros(N, 1);
u(F_nodes) = u_F;
u(B_nodes) = u_B;
u(D_nodes) = u_D;

%% STOP TIMER
t_FEM_BEM = toc;
disp(['FEM-BEM time: ' num2str(t_FEM_BEM) ' s'])

%% === DERIVATA NORMALE SU GAMMA (per post-processing e soluzione esterna) ===
% Segno più perché q_FEM = q_BEM
q_FEM = (G_bem \ H_bem) * u_B;   % ∂u/∂n con normale uscente da Omega_int

%% === DENSITA' DI CARICA SULLE PIASTRE (riga D) ===
% S_DF u_F + S_DD u_D + S_DB u_B = proiezione di ∂u/∂n sulle piastre
S_DF = S(D_nodes, F_nodes);
S_DD = S(D_nodes, D_nodes);
S_DB = S(D_nodes, B_nodes);

charge_projection = S_DF * u_F + S_DD * u_D + S_DB * u_B;
% charge_projection(i) = integral di (∂u/∂n * v_i) sulla piastra
% La densità di carica è sigma = e0 * ∂u/∂n

%% === POTENZIALE ESTERNO (formula di rappresentazione BEM) ===
% Per un punto y esterno a Gamma:
% u(y) = G_row(y) * q_BEM - H_row(y) * u_B
% dove q_BEM = -q_FEM
% Questa parte si implementa a parte valutando i kernel nel punto y

%% === POST-PROCESSING ===

% --- Campo elettrico per elemento ---
Ex = zeros(M,1); Ey = zeros(M,1);
for e = 1:M
    [Ex(e), Ey(e)] = Efield(e, x, y, conn, u);
end
Emag = sqrt(Ex.^2 + Ey.^2);


% --- Plot potenziale (mappa colori) ---
[up, xp, yp] = create2darray(x, y, u, delx, dely);

figure('Color', 'w');  % crea la figura con sfondo bianco
clf;                    % pulisce la figura
imagesc(xp, yp, up);
colormap(jet); axis equal tight;
set(gca,'Ydir','normal'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Potential (V)');
set(gcf,'Color',[1 1 1]); hold on
plot(x(ptop_node),    y(ptop_node),    'k', 'linewidth', 5)
plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 5)

% --- Plot potenziale (linee di livello) ---
figure('Color', 'w');  % crea la figura con sfondo bianco
clf;                    % pulisce la figura
contour(xp, yp, up, 'LevelStep', 0.1)
colormap(jet); axis equal tight; colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Potential (V) - contour');
set(gcf,'Color',[1 1 1]); hold on
plot(x(ptop_node),    y(ptop_node),    'k', 'linewidth', 5)
plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 5)

% --- Plot campo elettrico ---
[Emagp, xp, yp] = create2darray(xmid, ymid, Emag, delx, dely);

figure('Color', 'w');  % crea la figura con sfondo bianco
clf;                    % pulisce la figura
imagesc(xp, yp, Emagp);
colormap(jet); axis equal tight;
set(gca,'Ydir','normal'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Electric Field Intensity (V/m)');
set(gcf,'Color',[1 1 1]); hold on
quiver(xmid, ymid, Ex, Ey, 'color', 'k')
plot(x(ptop_node),    y(ptop_node),    'k', 'linewidth', 5)
plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 5)

% --- Plot derivata normale su Gamma ---
figure, clf
plot(q_FEM, 'b.-', 'MarkerSize', 10)
xlabel('Indice nodo su \Gamma');
ylabel('\partial u / \partial n');
title('Derivata normale su \Gamma (uscente da \Omega_{int})');
grid on


% Riordina q_FEM nell'ordine geometrico di B_ordered
q_ordered = q_FEM(perm);

figure('Color','w');
hold on
axis equal tight

% Plotta i segmenti di Gamma colorati per valore di q
n = length(B_ordered);
for k = 1:n
    n1 = k;
    n2 = mod(k, n) + 1;
    x1 = x(B_ordered(n1)); y1 = y(B_ordered(n1));
    x2 = x(B_ordered(n2)); y2 = y(B_ordered(n2));
    q_mid = 0.5*(q_ordered(n1) + q_ordered(n2));
    plot([x1 x2], [y1 y2], 'Color', ...
        interp1([-max(abs(q_ordered)) 0 max(abs(q_ordered))], ...
                [0 0 1; 1 1 1; 1 0 0], q_mid), ...
        'LineWidth', 4);
end

n = 256;
blue  = [0 0 1];
white = [1 1 1];
red   = [1 0 0];

cmap = [ ...
    linspace(blue(1),white(1),n/2)' ...
    linspace(blue(2),white(2),n/2)' ...
    linspace(blue(3),white(3),n/2)'; ...
    linspace(white(1),red(1),n/2)' ...
    linspace(white(2),red(2),n/2)' ...
    linspace(white(3),red(3),n/2)' ];

colormap(cmap)
cb = colorbar;
clim([-max(abs(q_ordered)) max(abs(q_ordered))]);
ylabel(cb, '\partial u / \partial n')
plot(x(ptop_node),    y(ptop_node),    'k', 'LineWidth', 4)
plot(x(pbottom_node), y(pbottom_node), 'k', 'LineWidth', 4)
xlabel('x (m)'); ylabel('y (m)')
title('\partial u / \partial n on \Gamma , E \cdot n , n positive going outside')

% --- Plot densità di carica sulle piastre ---
figure, clf
hold on
plot(x(ptop_node),    charge_projection(1:length(ptop_node)),             ...
     'r.-', 'MarkerSize', 10, 'DisplayName', 'Piastra superiore')
plot(x(pbottom_node), charge_projection(length(ptop_node)+1:end),         ...
     'b.-', 'MarkerSize', 10, 'DisplayName', 'Piastra inferiore')
xlabel('x (m)');
ylabel('$\int \frac{\partial u}{\partial n} \cdot v_i \, ds$', 'Interpreter', 'latex');
title('Proiezione densità di carica sulle piastre');
legend; grid on

% Visualizza distribuzione di E su d_elm
figure
histogram(Emag(d_elm), 50)
xlabel('|E| (V/m)'); ylabel('Conteggio elementi');
title('Distribuzione campo elettrico in d\_elm')

% Trova gli elementi problematici
disp("------- E Field -------")
E_teorico = (V1 - V2) / dy;
fattore_soglia = 5; % in percentuale
soglia = (1 + fattore_soglia/100) * E_teorico;
bad_elm  = d_elm(Emag(d_elm) > soglia);
good_elm = d_elm(Emag(d_elm) <= soglia);

disp(['E_teorico = ' num2str(E_teorico) ' V/m'])
disp(['Fattore Soglia = ' num2str(fattore_soglia) '%'])
disp(['E_soglia    = ' num2str(soglia)    ' V/m'])
disp(['Elementi con E > soglia: ' num2str(length(bad_elm))])
disp(['Elementi con E <= soglia: ' num2str(length(good_elm))])

% --- Flusso di E
fluxE = sum(q_FEM);
disp(['Flux of E on Gamma = ' num2str(fluxE) ' V*m'])


% Visualizza dove sono i bad elements

figure
triplot(conn, x, y, 'Color', [0.8 0.8 0.8]); hold on
if(length(bad_elm) > 0)
    triplot(conn(bad_elm,:), x, y, 'r');
end
plot(x(ptop_node), y(ptop_node), 'k', 'linewidth', 3)
plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 3)
title('Elementi con E anomalo (rosso)')


% --- Capacitanza ---
disp("------- Capacitanza -------")
% Totale (dominio FEM)
C = (1/(V1-V2)^2) * sum(eps_e .* (Ex.^2 + Ey.^2) .* area);
disp(['C (FEM region) = ' num2str(C*1e12) ' pF/m'])

% Interna al condensator
C_small = (1/(V1-V2)^2) * sum(eps_e(d_elm) .* (Ex(d_elm).^2 + Ey(d_elm).^2) .* area(d_elm));
disp(['C (inside) = ' num2str(C_small*1e12) ' pF/m'])

% Capacità escludendo elementi singolari
C_clean = (1/(V1-V2)^2) * sum(eps_e(good_elm) .* ...
          (Ex(good_elm).^2 + Ey(good_elm).^2) .* area(good_elm));
disp(['C (NO bad el) = ' num2str(C_clean*1e12) ' pF/m'])

% --- Carica Totale ---
disp("------- Carica Piastre -------")
n_top = length(ptop_node);
Q_top = sum(charge_projection(1:n_top));
Q_bot = sum(charge_projection(n_top+1:end));

disp(['Q piastra superiore = ' num2str(Q_top) ' C/m'])
disp(['Q piastra inferiore = ' num2str(Q_bot) ' C/m'])
disp(['Somma Q_top + Q_bot = ' num2str(Q_top + Q_bot) ' C/m  (deve essere ~0)'])
disp(['Rapporto Q_top/Q_bot = ' num2str(Q_top/Q_bot) '  (deve essere ~-1)'])


%% ERROR
error_FEM(dx,dy,delx,dely,e0,er,V1,V2,u,t_FEM_BEM);

%% SALVATAGGI
if(SALVATAGGIO == 1)

    % Fine del logging
    diary off


    % Salva tutte le figure
    %saveas(figure(1), fullfile(folder_name, 'mesh.png'));
    %saveas(figure(2), fullfile(folder_name, 'element_quality.png'));
    saveas(figure(1), fullfile(folder_name, 'potential_map.png'));
    saveas(figure(2), fullfile(folder_name, 'potential_contour.png'));
    saveas(figure(3), fullfile(folder_name, 'efield_map.png'));
    saveas(figure(4), fullfile(folder_name, 'q_gamma.png'));
    saveas(figure(5), fullfile(folder_name, 'q_gamma_ordered.png'));
    saveas(figure(6), fullfile(folder_name, 'charge_density.png'));
    saveas(figure(7), fullfile(folder_name, 'E_histogram.png'));
    saveas(figure(8), fullfile(folder_name, 'bad_elements.png'));
    saveas(figure(9), fullfile(folder_name, 'L2_norm.png'));
    saveas(figure(10), fullfile(folder_name, 'time.png'));
end




