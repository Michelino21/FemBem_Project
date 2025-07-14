function [errors,solutions,femregion,Dati]=C_main2D(TestName,nRef)
%==========================================================================
% Solution of the Poisson's problem with linear finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Matrices    : (struct) fe stiffness and mass matrices
%          Dati        : (struct)  see C_dati.m
%          
% Usage: 
%    [errors,solutions,femregion,Matrices,Dati] = C_main2D('Test1',3)
 


addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Dati = C_dati(TestName);
Dati.nRefinement = nRef;

%==========================================================================
% MESH GENERATION
%==========================================================================

[Region] = C_create_mesh(Dati);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = C_create_femregion(Dati,Region); 

%==========================================================================
% BUILD FINITE ELEMENT MATRICES (M and A) and RIGHT-HAND SIDE (b)
%==========================================================================

[M_no_bc,A_no_bc,r_no_bc] = C_matrix2D(Dati,femregion);

% time integration parameters
t = 0;              %initial time
T = Dati.T;         %final time
dt = Dati.dt;       %time-step
theta = Dati.theta; %theta-method

x = femregion.coord(:,1);
y = femregion.coord(:,2);

% Evaluation of initial condition
u0 = eval(Dati.initialcond);

% Snapshot of the initial condition
if Dati.visual_graph == 'Y'
    C_snapshot_sol(femregion, u0, Dati, t)
end

% evaluation of force at time tk
f0 = eval(Dati.forcetime);
K_no_bc = M_no_bc+dt*theta*A_no_bc;
W = M_no_bc-dt*(1-theta)*A_no_bc;
for t = dt : dt : T 
    f1 = eval(Dati.forcetime);    
    rhs_no_bc = W*u0 + dt*(theta*f1+(1-theta)*f0)*r_no_bc;
    [K,rhs,u_g] = C_bound_cond2D(K_no_bc,rhs_no_bc,femregion,Dati,t);
    u1 = K \ rhs + u_g;
    u0 = u1;
    f0 = f1;
    if Dati.visual_graph == 'Y' && mod(round(t/dt), 10) == 0
        figure()
        C_snapshot_sol(femregion, u0, Dati, t)
    end
end

%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions] = C_postprocessing(Dati,femregion,u1,t);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
errors = [];
if (Dati.plot_errors)
    [errors] = C_compute_errors(Dati,femregion,solutions,t);
end



