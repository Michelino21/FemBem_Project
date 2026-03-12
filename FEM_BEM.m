%% FEM_BEM.m
% Coupling FEM-BEM per condensatore a facce piane parallele 2D
clear all; close all; clc;
addpath(genpath('src'));

SALVATAGGIO = 0';

%% SALVATAGGI
if(SALVATAGGIO == 1)

    % Scegli il nome della cartella
    folder_name = './Documentation/figures/coarse_2';  % cambia qui per ogni run
    
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
Lx   = 0.10;          % lunghezza piastre
Ly   = 0.01;          % separazione piastre  →  rapporto Lx/Ly = 10
delx = 0.005;         % passo mesh
dely = delx;

margine = 3 * Ly;     % = 0.03 m per lato
X_mesh = Lx + 2*margine;  % = 0.16 m
Y_mesh = Ly + 2*margine;  % = 0.07 m


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
    meshgen_cap(X_mesh, Y_mesh, Lx, Ly, delx, dely);    
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

%%Risoluzioen (G^-1)*H
BEM_matrix = G_bem \ H_bem;

% Nota: segno + perché n_FEM = n_BEM
K = [S_FF,                        S_FB; ...
     S_BF,   S_BB - T_BB * (BEM_matrix)];

rhs = -[S_FD; S_BD] * u_D;

%% === SOLUZIONE ===
sol   = K \ rhs;
u_F   = sol(1:length(F_nodes));
u_B   = sol(length(F_nodes)+1:end);

% Ricostruzione del potenziale globale (Step 1)
u = zeros(N, 1);
u(F_nodes) = u_F;
u(B_nodes) = u_B;
u(D_nodes) = u_D;

%% STOP TIMER
t_FEM_BEM = toc;
disp(['FEM-BEM time: ' num2str(t_FEM_BEM) ' s'])

%% ====== POST-PROCESSING ===================
%
%============================================
display("POST-PROCESSING");
% Step 1.a: Potenziale Interno
% Step 1.b: Campo Elettrico Interno
% Step 2: Derivata Normale su Gamma
% Step 3: Potenziale Esterno
% Step 4.a: Densita superficiale di carica
% Step 4.b: Capacità
% Step 5: Analisi elementi di campo elettrico problematici.

%% === POTENZIALE INTERNO === (Step 1.a)

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


%% === CAMPO ELETTRICO INTERNO === (Step 1.b)

% --- Campo elettrico per elemento ---
Ex = zeros(M,1); Ey = zeros(M,1);
for e = 1:M
    [Ex(e), Ey(e)] = Efield(e, x, y, conn, u);
end
Emag = sqrt(Ex.^2 + Ey.^2);

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


% Visualizza distribuzione di E su d_elm
figure
histogram(Emag(d_elm), 50)
xlabel('|E| (V/m)'); ylabel('Conteggio elementi');
title('Distribuzione campo elettrico in d\_elm')

% Trova gli elementi problematici
disp("------- E Field -------")
E_teorico = (V1 - V2) / Ly;
fattore_soglia = 5; % in percentuale
soglia = (1 + fattore_soglia/100) * E_teorico;
bad_elm  = d_elm(Emag(d_elm) > soglia);
good_elm = d_elm(Emag(d_elm) <= soglia);

disp(['E_teorico = ' num2str(E_teorico) ' V/m'])
disp(['Fattore Soglia = ' num2str(fattore_soglia) '%'])
disp(['E_soglia    = ' num2str(soglia)    ' V/m'])
disp(['Elementi con E > soglia: ' num2str(length(bad_elm))])
disp(['Elementi con E <= soglia: ' num2str(length(good_elm))])


% Visualizza dove sono i bad elements

figure
triplot(conn, x, y, 'Color', [0.8 0.8 0.8]); hold on
if(length(bad_elm) > 0)
    triplot(conn(bad_elm,:), x, y, 'r');
end
plot(x(ptop_node), y(ptop_node), 'k', 'linewidth', 3)
plot(x(pbottom_node), y(pbottom_node), 'k', 'linewidth', 3)
title('Elementi con E anomalo (rosso)')


%% === DERIVATA NORMALE SU GAMMA (per post-processing e soluzione esterna) === (Step 2)
% Segno più perché q_FEM = q_BEM
q_FEM = (BEM_matrix) * u_B;   % ∂u/∂n con normale uscente da Omega_int


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


% --- Flusso di E
fluxE = sum(q_FEM);
disp(['Flux of E on Gamma = ' num2str(fluxE) ' V*m'])

%% === POTENZIALE ESTERNO (formula di rappresentazione BEM) === (Step 3)
% Per un punto y esterno a Gamma:
% u(y) = G_row(y) * q_BEM - H_row(y) * u_B
% dove q_BEM = -q_FEM
% Questa parte si implementa a parte valutando i kernel nel punto y

display("====== POTENZIALE ESTERNO ======");

% ========================================================================
% STEP 2: Definisci il centro di riferimento del condensatore
% ========================================================================
x_plate_min = -Lx/2;
x_plate_max = Lx/2;
y_plate_top = Ly/2;
y_plate_bottom = -Ly/2;

x_center = 0;
y_center = 0;

% Dimensioni caratteristiche
char_size = max(Lx, Ly);

%fprintf('Centro condensatore: (%.4f, %.4f)\n', x_center, y_center);
%fprintf('Dimensione caratteristica: %.4f m\n', char_size);

% ========================================================================
% STEP 3: Funzioni kernel per BEM
% ========================================================================
% Kernel fondamentale G per Laplace 2D
G_kernel = @(x1, y1, x2, y2) -(1/(2*pi)) * log(sqrt((x2-x1).^2 + (y2-y1).^2) + 1e-16);

% Derivata normale del kernel H
H_kernel = @(x1, y1, x2, y2, nx, ny) -(1/(2*pi)) * ...
    ((x1-x2).*nx + (y1-y2).*ny) ./ ((x1-x2).^2 + (y1-y2).^2 + 1e-16);

% ========================================================================
% STEP 4: Calcola normali ai nodi di Γ (se non già disponibili)
% ========================================================================
% Usa i nodi ordinati per calcolare le normali
N_B = length(B_ordered);
normal_x = zeros(N_B, 1);
normal_y = zeros(N_B, 1);

for i = 1:N_B
    % Nodi precedente e successivo (percorso chiuso)
    i_prev = mod(i - 2, N_B) + 1;
    i_next = mod(i, N_B) + 1;
    
    node_prev = B_ordered(i_prev);
    node_curr = B_ordered(i);
    node_next = B_ordered(i_next);
    
    % Vettore tangente (media tra i due segmenti)
    tx1 = x(node_curr) - x(node_prev);
    ty1 = y(node_curr) - y(node_prev);
    tx2 = x(node_next) - x(node_curr);
    ty2 = y(node_next) - y(node_curr);
    
    tx = tx1 + tx2;
    ty = ty1 + ty2;
    t_norm = sqrt(tx^2 + ty^2);
    
    if t_norm > 1e-10
        tx = tx / t_norm;
        ty = ty / t_norm;
    end
    
    % Normale (ruotata 90° in senso orario per normale esterna)
    normal_x(i) = ty;
    normal_y(i) = -tx;
end

% Riporta le normali all'ordine di B_nodes
[~, perm_back] = ismember(B_nodes, B_ordered);
nx_B = normal_x(perm_back);
ny_B = normal_y(perm_back);

% ========================================================================
% STEP 5: Definisci punti esterni lungo 4 direzioni radiali
% ========================================================================
n_points = 50;  % Numero di punti per direzione
r_min = 2.0 * char_size;   % Distanza minima (fuori dal dominio FEM)
r_max = 10.0 * char_size;  % Distanza massima
r_values = linspace(r_min, r_max, n_points);

% 4 direzioni: +y, -y, +x, -x
directions = {
    [0, 1],   'Sopra (+y)',    'b';  % Direzione, Label, Colore
    [0, -1],  'Sotto (-y)',    'r';
    [1, 0],   'Destra (+x)',   'g';
    [-1, 0],  'Sinistra (-x)', 'm'
};

n_dir = size(directions, 1);

% Matrice per salvare i risultati
u_ext = zeros(n_points, n_dir);

% ========================================================================
% STEP 6: Calcola u esterno usando BEM representation formula
% ========================================================================
fprintf('Calcolo potenziale esterno lungo %d direzioni...\n', n_dir);

for d = 1:n_dir
    dir_vec = directions{d, 1};
    
    for i = 1:n_points
        r = r_values(i);
        
        % Punto esterno
        y_ext = x_center + r * dir_vec(1);
        x_ext = y_center + r * dir_vec(2);  % Nota: y qui è la coord verticale
        
        % Calcola i vettori G(y_ext, x_B) e H(y_ext, x_B)
        G_vec = zeros(length(B_nodes), 1);
        H_vec = zeros(length(B_nodes), 1);
        
        for j = 1:length(B_nodes)
            node_j = B_nodes(j);
            
            G_vec(j) = G_kernel(x_ext, y_ext, y(node_j), x(node_j));
            H_vec(j) = H_kernel(x_ext, y_ext, y(node_j), x(node_j), ...
                                ny_B(j), nx_B(j));
        end
        
        % Formula BEM: u(y) = G·q - H·u_B
        u_ext(i, d) = G_vec' * q_FEM - H_vec' * u_B;
    end
    
    fprintf('  Direzione %s: completata\n', directions{d, 2});
end

% ========================================================================
% STEP 7: Fit con legge 1/r (per dipolo 2D)
% ========================================================================
fit_coeff = zeros(n_dir, 1);

for d = 1:n_dir
    % Fit lineare in log-log: log(u) = log(A) - log(r)
    % Usa solo punti lontani (r > 5*char_size) per evitare effetti near-field
    idx_far = r_values > 5 * char_size;
    
    if sum(idx_far) > 5
        p = polyfit(log(r_values(idx_far)), log(abs(u_ext(idx_far, d))), 1);
        fit_coeff(d) = exp(p(2));  % A = exp(intercetta)
        slope = p(1);
        fprintf('  %s: pendenza = %.2f (teorica: -1 per dipolo)\n', ...
                directions{d, 2}, slope);
    end
end

% ========================================================================
% GRAFICO 1: Plot Radiali log-log
% ========================================================================
figure('Position', [100, 100, 1200, 800]);

for d = 1:n_dir
    subplot(2, 2, d);
    
    % Plot dei dati BEM
    loglog(r_values, abs(u_ext(:, d)), 'o', ...
           'Color', directions{d, 3}, 'MarkerSize', 6, 'LineWidth', 1.5);
    hold on;
    
    % Fit 1/r
    u_fit = fit_coeff(d) ./ r_values;
    loglog(r_values, u_fit, '--', 'Color', directions{d, 3}, 'LineWidth', 2);
    
    % Linea di riferimento pendenza -1
    r_ref = [r_min, r_max];
    u_ref = u_ext(1, d) * (r_min ./ r_ref);
    loglog(r_ref, abs(u_ref), 'k:', 'LineWidth', 1);
    
    grid on;
    xlabel('Distanza r [m]');
    ylabel('|u(r)| [V]');
    title(sprintf('Direzione: %s', directions{d, 2}));
    legend('BEM', 'Fit A/r', 'Pendenza -1', 'Location', 'best');
    hold off;
end

sgtitle('Potenziale Esterno vs Distanza (Log-Log)');

% ========================================================================
% GRAFICO 2: Mappa 2D del potenziale esterno
% ========================================================================
figure('Position', [100, 100, 900, 700]);

% Griglia di punti esterni
n_grid = 100;
x_grid_min = x_center - 8 * char_size;
x_grid_max = x_center + 8 * char_size;
y_grid_min = y_center - 8 * char_size;
y_grid_max = y_center + 8 * char_size;

[X_grid, Y_grid] = meshgrid(linspace(x_grid_min, x_grid_max, n_grid), ...
                             linspace(y_grid_min, y_grid_max, n_grid));

U_grid = zeros(size(X_grid));

fprintf('Calcolo mappa 2D (%d × %d punti)...\n', n_grid, n_grid);

for i = 1:n_grid
    for j = 1:n_grid
        y_ext = X_grid(i, j);
        x_ext = Y_grid(i, j);
        
        % Salta punti troppo vicini al boundary (dentro il dominio FEM)
        r_min_bound = 1.5 * char_size;
        if sqrt((y_ext - x_center)^2 + (x_ext - y_center)^2) < r_min_bound
            U_grid(i, j) = NaN;
            continue;
        end
        
        % Calcola BEM
        G_vec = zeros(length(B_nodes), 1);
        H_vec = zeros(length(B_nodes), 1);
        
        for k = 1:length(B_nodes)
            node_k = B_nodes(k);
            G_vec(k) = G_kernel(x_ext, y_ext, y(node_k), x(node_k));
            H_vec(k) = H_kernel(x_ext, y_ext, y(node_k), x(node_k), ...
                                ny_B(k), nx_B(k));
        end
        
        U_grid(i, j) = G_vec' * q_FEM - H_vec' * u_B;
    end
    if mod(i, 10) == 0
        fprintf('  Progresso: %d/%d righe\n', i, n_grid);
    end
end

% Plot contour
contourf(X_grid, Y_grid, U_grid, 20, 'LineWidth', 0.5);
hold on;
colorbar;
colormap(jet);

% Aggiungi il boundary Γ
plot(x(B_ordered), y(B_ordered), 'k-', 'LineWidth', 2);

% Aggiungi le piastre
plot(x(ptop_node), y(ptop_node), 'w-', 'LineWidth', 3);
plot(x(pbottom_node), y(pbottom_node), 'w-', 'LineWidth', 3);

% Cerchi di riferimento
theta = linspace(0, 2*pi, 100);
for r_circle = [3, 5, 7] * char_size
    plot(x_center + r_circle * cos(theta), ...
         y_center + r_circle * sin(theta), 'w--', 'LineWidth', 1);
end

axis equal;
xlabel('x [m]');
ylabel('y [m]');
title('Mappa 2D: Potenziale Esterno u(x,y)');
hold off;

disp('Grafici completati!');



%% === DENSITA' DI CARICA SULLE PIASTRE  === (Step 4)

% ========================================================================
% Recupera E su piastra alta e bassa
% ========================================================================


% Trova elementi che contengono almeno un nodo della piastra superiore
elements_top = [];
for e = 1:M
    nodes_in_element = conn(e, :);
    if any(ismember(nodes_in_element, ptop_node))
        % Verifica che l'elemento sia nel gap (sotto la piastra)
        y_centroid = mean(y(nodes_in_element));
        y_plate_top = y(ptop_node(1));  % y della piastra superiore
        if y_centroid < y_plate_top  % Elemento sotto la piastra
            elements_top = [elements_top; e];
        end
    end
end

% Trova elementi adiacenti alla piastra inferiore
elements_bottom = [];
for e = 1:M
    nodes_in_element = conn(e, :);
    if any(ismember(nodes_in_element, pbottom_node))
        % Verifica che l'elemento sia nel gap (sopra la piastra)
        y_centroid = mean(y(nodes_in_element));
        y_plate_bottom = y(pbottom_node(1));  % y della piastra inferiore
        if y_centroid > y_plate_bottom  % Elemento sopra la piastra
            elements_bottom = [elements_bottom; e];
        end
    end
end

% ========================================================================
% Fai la proiezione sulla normale
% ========================================================================

% Normale alla piastra superiore (verso il basso, verso il gap)
n_top = [0; -1];

% Normale alla piastra inferiore (verso l'alto, verso il gap)
n_bottom = [0; 1];

% Proiezione del campo E sulla normale per piastra superiore
En_top = zeros(length(elements_top), 1);
for i = 1:length(elements_top)
    e = elements_top(i);
    E_vec = [Ex(e); Ey(e)];
    En_top(i) = dot(E_vec, n_top);
end

% Proiezione del campo E sulla normale per piastra inferiore
En_bottom = zeros(length(elements_bottom), 1);
for i = 1:length(elements_bottom)
    e = elements_bottom(i);
    E_vec = [Ex(e); Ey(e)];
    En_bottom(i) = dot(E_vec, n_bottom);
end

% ========================================================================
% Calcola la densità di carica per elemento
% ========================================================================

% Densità di carica per elementi adiacenti alla piastra superiore
sigma_top = e0 * En_top;

% Densità di carica per elementi adiacenti alla piastra inferiore
sigma_bottom = e0 * En_bottom;

% ========================================================================
% Calcola Q (carica totale per unità di lunghezza)
% ========================================================================

% Per ogni elemento, calcola l'area (o lunghezza del lato sulla piastra)
% Assumendo elementi triangolari in 2D

% Carica sulla piastra superiore
Q_top = 0;
for i = 1:length(elements_top)
    e = elements_top(i);
    nodes_in_element = conn(e, :);
    
    % Trova il lato dell'elemento che giace sulla piastra
    % (i nodi che appartengono a ptop_node)
    nodes_on_plate = nodes_in_element(ismember(nodes_in_element, ptop_node));
    
    if length(nodes_on_plate) >= 2
        % Calcola lunghezza del lato sulla piastra
        x1 = x(nodes_on_plate(1));
        y1 = y(nodes_on_plate(1));
        x2 = x(nodes_on_plate(2));
        y2 = y(nodes_on_plate(2));
        
        length_segment = sqrt((x2-x1)^2 + (y2-y1)^2);
        
        % Contributo alla carica totale
        Q_top = Q_top + sigma_top(i) * length_segment;
    end
end

% Carica sulla piastra inferiore
Q_bottom = 0;
for i = 1:length(elements_bottom)
    e = elements_bottom(i);
    nodes_in_element = conn(e, :);
    
    % Trova il lato dell'elemento che giace sulla piastra
    nodes_on_plate = nodes_in_element(ismember(nodes_in_element, pbottom_node));
    
    if length(nodes_on_plate) >= 2
        % Calcola lunghezza del lato sulla piastra
        x1 = x(nodes_on_plate(1));
        y1 = y(nodes_on_plate(1));
        x2 = x(nodes_on_plate(2));
        y2 = y(nodes_on_plate(2));
        
        length_segment = sqrt((x2-x1)^2 + (y2-y1)^2);
        
        % Contributo alla carica totale
        Q_bottom = Q_bottom + sigma_bottom(i) * length_segment;
    end
end

% ========================================================================
% Risultati
% ========================================================================

fprintf('Carica totale piastra superiore: Q_top = %.6e C/m\n', Q_top);
fprintf('Carica totale piastra inferiore: Q_bottom = %.6e C/m\n', Q_bottom);
fprintf('Somma (dovrebbe essere ~0): Q_tot = %.6e C/m\n', Q_top + Q_bottom);
fprintf(['Rapporto Q_top/Q_bot = ' num2str(Q_top/Q_bottom) '  (deve essere ~-1)'])

% Densità di carica media
sigma_top_mean = mean(sigma_top);
sigma_bottom_mean = mean(sigma_bottom);

fprintf('Densità di carica media piastra superiore: %.6e C/m^2\n', sigma_top_mean);
fprintf('Densità di carica media piastra inferiore: %.6e C/m^2\n', sigma_bottom_mean);

% Opzionale: visualizza la distribuzione di carica
figure;
subplot(2,1,1);
plot(sigma_top, 'o-');
title('Densità di carica - Piastra Superiore');
xlabel('Indice elemento');
ylabel('\sigma [C/m^2]');
grid on;

subplot(2,1,2);
plot(sigma_bottom, 'o-');
title('Densità di carica - Piastra Inferiore');
xlabel('Indice elemento');
ylabel('\sigma [C/m^2]');
grid on;


%%  === CAPACITA'  === (Step 4.b)

% Differenza di potenziale
Delta_V = V1 - V2;

% Capacitanza usando la carica dalla piastra superiore
C_from_top = abs(Q_top) / Delta_V;

% Capacitanza usando la carica dalla piastra inferiore
C_from_bottom = abs(Q_bottom) / Delta_V;

% Capacitanza media (dovrebbero essere quasi uguali)
C = (C_from_top + C_from_bottom) / 2;

% ========================================================================
% Confronto con soluzione analitica (condensatore piano ideale)
% ========================================================================

% Capacitanza teorica per condensatore piano infinito
C_theoretical = e0 * er * Lx/Ly;

% Errore percentuale
error_percent = abs(C - C_theoretical) / C_theoretical * 100;


fprintf('\n========== CAPACITANZA ==========\n');
fprintf('Capacitanza (da piastra superiore): C = %.6e F/m\n', C_from_top);
fprintf('Capacitanza (da piastra inferiore): C = %.6e F/m\n', C_from_bottom);
fprintf('Capacitanza media: C = %.6e F/m\n', C);
fprintf('Capacitanza teorica (piano infinito): C_th = %.6e F/m\n', C_theoretical);
fprintf('Errore relativo: %.2f%%\n', error_percent);
fprintf(' ');
fprintf('Lx = %.6f m\n', Lx);
fprintf('Ly = %.6f m\n', Ly);
fprintf('Q_{top}    = %+.6e C/m\n', Q_top);
fprintf('Q_{bottom} = %+.6e C/m\n', Q_bottom);
fprintf('Q_{sum} = %+.6e C/m\n', Q_top + Q_bottom);
fprintf('ΔV         = %.6f V\n', Delta_V);


% USATO PER QUANDO CALCOLAVO CAPACITANZA CON ENERGIA.
% disp("------- Capacitanza -------")
% % Totale (dominio FEM)
% C = (1/(V1-V2)^2) * sum(eps_e .* (Ex.^2 + Ey.^2) .* area);
% disp(['C (FEM region) = ' num2str(C*1e12) ' pF/m'])
% 
% % Interna al condensator
% C_small = (1/(V1-V2)^2) * sum(eps_e(d_elm) .* (Ex(d_elm).^2 + Ey(d_elm).^2) .* area(d_elm));
% disp(['C (inside) = ' num2str(C_small*1e12) ' pF/m'])
% 
% % Capacità escludendo elementi singolari
% C_clean = (1/(V1-V2)^2) * sum(eps_e(good_elm) .* ...
%           (Ex(good_elm).^2 + Ey(good_elm).^2) .* area(good_elm));
% disp(['C (NO bad el) = ' num2str(C_clean*1e12) ' pF/m'])



%% ERROR
error_FEM(Lx,Ly,delx,dely,e0,er,V1,V2,u,t_FEM_BEM);

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




disp('FINE CODICE');