%% dO2

% bio
Jo2_tmp = ro2p*(-prod_p) + ro2p*(remin_p + regen_p + tr.dop/taup).*(tr.O2./(tr.O2 + KO2));

% check for low o2
io2_islo = find( Jo2_tmp<0 & tr.O2 <= (o2_crit + eps_lim));
io2_golo = find( Jo2_tmp<0 & tr.O2 > (o2_crit + eps_lim) & (-Jo2_tmp)*dt > (tr.O2 - o2_crit - eps_lim) );
% correct bio term
Jo2 = Jo2_tmp;
Jo2(io2_islo) = 0;
Jo2(io2_golo) = - (tr.O2(io2_golo) - o2_crit - eps_lim)/dt;
res_oxid = Jo2_tmp - Jo2;
res_oxid_p = res_oxid/ro2p;
%res_oxid_p(isub) = 0;
%res_oxid_p(grid.ZT3d(iocn) < 150) = 0;

% gas exchange
Fo2 = obs.kw_o2.*( obs.o2sat - tr.O2 );
gxo2 = Fo2./grid.dzt(1);