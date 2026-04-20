%% dFe


if ~do.POM;
    exp_fe = ef.*mort_fe;
    remin_fe = Qrem*(exp_fe);
    regen_fe = (1 - ef).*mort_fe;

else % do.POM
    remin_fe = tr.pofe/taupom;
    regen_fe = (1 - phie)*mort_fe;
    Jpofe = Qsink*tr.pofe + phie*mort_fe - tr.pofe/taupom;
end

Jfe = -prod_fe + remin_fe + regen_fe;

% deposition
Sfe_dep = do.Fe_sms*1000*(alpha_fe*fe_dep/56)./grid.DZT3d(iocn); % mmol/m3/yr

% find free fe
[fe_free,junk1,junk2] = solve_Fe(tr.Fe,Ligand,beta_fe); % mmol/m3

% scavenging
Sfe_scav = do.Fe_sms*ksc*fe_free;  % /yr
%Sfe_scav(isub) = 0;
