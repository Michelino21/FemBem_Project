function [u_ex]=C_eval_exact_sol(femregion, exact_sol,t)
%% [u_ex]=C_eval_exact_sol(femregion, exact_sol,t)
%==========================================================================
% EVALUATION OF THE EXACT SOLUTION ON THE GRID POINTS
%==========================================================================
%    called in C_postprocessing.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          exact_sol   : (string)  expression of exact solution
%          t           : (real) current time
%
%    OUTPUT:
%          u_ex        : (spare(ndof,1) real) exact solution vector
%



dof = femregion.dof;
x = dof(:,1);
y = dof(:,2);
u_ex = eval(exact_sol);
    

