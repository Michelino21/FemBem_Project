function C_snapshot_sol(femregion, uh, Dati, time)
%% C_snapshot_sol(femregion, uh, Dati, time)
%==========================================================================
% PLOT THE EXACT SOLUTION ON THE DOFS
%==========================================================================
%    called in C_matrix2D.m
%
%    INPUT:
%          femregion   : (struct)  see C_create_femregion.m
%          uh          : (sparse(ndof,1) real) solution vector
%          Dati        : (struct) see C_Dati.m
%          time        : (real) current time


dof=femregion.dof;

x1=femregion.domain(1,1);
x2=femregion.domain(1,2);
y1=femregion.domain(2,1);
y2=femregion.domain(2,2);

% M=max(uh);
% m=min(uh);
M=1;
m=-1;
% if (abs(m-M) < 0.1)
%     M=m+1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VISUALIZZAZIONE DELLE SOLUZIONI MESHATE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
trisurf(femregion.connectivity(1:3,:)',femregion.coord(:,1),femregion.coord(:,2),full(uh));
title(['u_h(x,y) at time: ' num2str(time)]);
xlabel('x-axis'); 
ylabel('y-axis');
axis([x1,x2,y1,y2,m,M]); cc = colorbar;
pause(0.1);




