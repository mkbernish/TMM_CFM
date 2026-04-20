% tracers
tr_str = 'tr'; ntr = 0;
if do.P; tr_str = [tr_str '_P']; ntr = ntr+1; end 
if do.N; tr_str = [tr_str '_N']; ntr = ntr+1; end
if do.Fe; tr_str = [tr_str '_Fe']; ntr = ntr+1; end
if do.O2; tr_str = [tr_str '_O2']; ntr = ntr+1; end
if do.phyto; tr_str = [tr_str '_phy']; ntr = ntr+1; end
if do.diaz; tr_str = [tr_str '_diaz']; ntr = ntr+1; end
if ntr==6; tr_str = 'tr_all'; end
if do.DOM; tr_str = [tr_str '_DOM']; end
if do.cons_fe
    fe_str = '_consFe';
else
    fe_str = 'varFe';
end

% denit
if do.denit; 
    if do.fixed_denit; D_str = ['_Dfix_' num2str(Dtot)]; 
    else D_str = ['_Dfree']; end
else D_str = [];
end

% diazotrophs
if do.diaz & do.Fe;
    diazFe_str = ['_diazFe_' num2str(Fe_fac)];
else diazFe_str = [];
end

% ro 
if do.var_n2p
    if do.spec_rost
        ro_str = ['_rst' num2str(ro_st)];
    else
        ro_str = ['_rdi' num2str(ro_di)];
    end
else
    ro_str = [];
end

% glacial
if do.glacial_O2 | do.glacial_Fe | do.glacial_SED;
    G_str = '_Glacial';
    
    if do.glacial_O2==1;
        G_str = [G_str '_O2'];
    elseif do.glacial_O2==5;
        G_str = [G_str '_5O2']; 
    else G_str = G_str; end

    if do.glacial_Fe;
        G_str = [G_str '_Fe']; end
    
    if do.glacial_SED==1;
        G_str = [G_str '_SEDw'];
    elseif do.glacial_SED==2;
        G_str = [G_str '_SEDs'];
    else G_str = G_str; end
    
else
    G_str = [];
end

% phyto
if do.phyto & ~do.qssa
    phy_str = '_explic';
else
    phy_str = [];
end
if fix_cost~=0.85;
    phy_str = [phy_str '_fc' num2str(fix_cost*100)];
end

% O2
if do.O2;
    if KO2==0;
        o2_str=[];
    else
        o2_str=['_KO2_' num2str(KO2)];
    end
end
if do.coupled_fdi
    fdi_str = ['_cfdi_'];
else
    fdi_str = ['_ofdi_'];
end

if do.coupled_mub
    mub_str = ['_cmub_'];
else
    mub_str = ['_omub_'];
end

% full name
nmod = [tr_str D_str ro_str diazFe_str phy_str G_str fdi_str mub_str fe_str];