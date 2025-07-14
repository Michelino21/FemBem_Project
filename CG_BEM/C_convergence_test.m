TestName = 'Test1';
Dati = C_dati(TestName);

nrefs = Dati.refinement_vector;

for i = 1:length(nrefs)
    [errors,solutions,femregion,Dati]=C_main2D(TestName,nrefs(i));
    h(i) = femregion.h;
    err_L2(i) = errors.Error_L2;
    err_H1(i) = errors.Error_H1;
end

%% 

err_L2
err_H1