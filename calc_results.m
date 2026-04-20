%% Calculate Results %%

%% PO4 & DOP
if do.P
    Pmat = 0*M3d+NaN; Pmat(iocn) = tr.P;
    Jpmat = 0*M3d+NaN; Jpmat(iocn) = Jp;

    if do.POM;
        POPmat = M3d*0 + NaN; POPmat(iocn) = tr.pop;
        POPsink = M3d*0 + NaN; POPsink(iocn) = tr.pop.*wsink;
        expP = POPsink(:,:,kprod);
        eP = nansum(nansum(expP.*AREA(:,:,kprod)))/1000;
    else
        expP3d = 0*M3d+NaN; expP3d(iocn) = exp_p;
        expP = nansum(expP3d(:,:,1:kprod).*grid.DZT3d(:,:,1:kprod),3); expP(M3d(:,:,1)==0)=NaN;
        eP = (exp_p)'*dv/1000; % in molP/yr
    end
    eC = eP*106*12/1e15; % in GTC/yr

    if do.DOM;
        DOPmat = 0*M3d+NaN; DOPmat(iocn) = tr.dop; end

end

%% NO3 & DON
if do.N
    Nmat = 0*M3d+NaN; Nmat(iocn) = tr.N;
    Jnmat = 0*M3d+NaN; Jnmat(iocn) = Jn;

    if do.POM;
        PONmat = M3d*0 + NaN; PONmat(iocn) = tr.pon;
        PONsink = M3d*0 + NaN; PONsink(iocn) = tr.pon.*wsink;
        expN = PONsink(:,:,kprod);
        eN = nansum(nansum(expN.*AREA(:,:,kprod)));
    else
        expN3d = 0*M3d+NaN; expN3d(iocn) = exp_n;
        expN = nansum(expN3d(:,:,1:kprod).*grid.DZT3d(:,:,1:kprod),3); expN(M3d(:,:,1)==0)=NaN;
        eN = (exp_n)'*dv/1000; % in molP/yr
    end

    if do.diaz;
        JNfixmat = 0*M3d+NaN; JNfixmat(iocn) = prod_f*rf;
        Nfix_2d = sum(JNfixmat.*grid.DZT3d,3,'omitnan');
        Nfix_2d(M3d(:,:,1)==0)=NaN;

        % basin integrals
        Nfix_ATL = (Nfix_2d(isurf).*MSKS.ATL(isurf))'*da(ivsurf)*14/1e15;
        Nfix_PAC = (Nfix_2d(isurf).*MSKS.PAC(isurf))'*da(ivsurf)*14/1e15;
        Nfix_IND = (Nfix_2d(isurf).*MSKS.IND(isurf))'*da(ivsurf)*14/1e15;
    end

    if do.DOM;
        DONmat = 0*M3d+NaN; DONmat(iocn) = tr.don; end

    if do.denit
        WCmat = 0*M3d+NaN; WCmat(iocn) = S_wc;
        SEDmat = 0*M3d+NaN; SEDmat(iocn) = S_sed;
        WC_2d = sum(WCmat.*grid.DZT3d,3,'omitnan'); WC_2d(M3d(:,:,1)==0)=NaN;
        SED_2d = sum(SEDmat.*grid.DZT3d,3,'omitnan'); SED_2d(M3d(:,:,1)==0)=NaN;
        WC_prof = squeeze(sum(sum(WCmat.*VOL,1,'omitnan'),2,'omitnan'))*14/1e15;
        SED_prof = squeeze(sum(sum(SEDmat.*VOL,1,'omitnan'),2,'omitnan'))*14/1e15;
        F_wc = sink_WC./(sink_WC + sink_SED);

        % basin integrals
        Den_ATL = ( (WC_2d(isurf) + SED_2d(isurf)).*MSKS.ATL(isurf) )'*da(ivsurf)*14/1e15;
        Den_PAC = ( (WC_2d(isurf) + SED_2d(isurf)).*MSKS.PAC(isurf) )'*da(ivsurf)*14/1e15;
        Den_IND = ( (WC_2d(isurf) + SED_2d(isurf)).*MSKS.IND(isurf) )'*da(ivsurf)*14/1e15;
    end

    % P* and N*
    Ps_obs = po4obs - no3obs/16;
    Ps = Pmat - Nmat/16;

    % create low latitude mask
    LMASK = M3d*0;
    LMASK(grid.YT3d<=35 & grid.YT3d>=-35) = 1;

    % basin mean P*
    Ps_ATL = ((MSKS.ATL(isurf).*LMASK(isurf).*Ps(isurf))'*dv(ivsurf))/((MSKS.ATL(isurf).*LMASK(isurf))'*dv(ivsurf));
    Ps_PAC = ((MSKS.PAC(isurf).*LMASK(isurf).*Ps(isurf))'*dv(ivsurf))/((MSKS.PAC(isurf).*LMASK(isurf))'*dv(ivsurf));
    Ps_obs_ATL = ((MSKS.ATL(isurf).*LMASK(isurf).*Ps_obs(isurf))'*dv(ivsurf))/((MSKS.ATL(isurf).*LMASK(isurf))'*dv(ivsurf));
    Ps_obs_PAC = ((MSKS.PAC(isurf).*LMASK(isurf).*Ps_obs(isurf))'*dv(ivsurf))/((MSKS.PAC(isurf).*LMASK(isurf))'*dv(ivsurf));

    % Pacific - Atlantic gradient
    Ps_grad = Ps_PAC - Ps_ATL;
    Ps_obs_grad = Ps_obs_PAC - Ps_obs_ATL;
end

%% Phytoplankton
if do.phyto
    Bmat = 0*M3d+NaN; Bmat(iocn) = tr.B;
    Fmat = 0*M3d+NaN; Fmat(iocn) = tr.F;
    Fixmat = 0*M3d+NaN; Fixmat(iocn) = S_nfix;
    QnSmat = 0*M3d+NaN; QnSmat(iocn) = tr.QnS;
    QpSmat = 0*M3d+NaN; QpSmat(iocn) = tr.QpS;
    QfeSmat = 0*M3d+NaN; QfeSmat(iocn) = tr.QfeS;
    QnLmat = 0*M3d+NaN; QnLmat(iocn) = tr.QnL;
    QpLmat = 0*M3d+NaN; QpLmat(iocn) = tr.QpL;
    QfeLmat = 0*M3d+NaN; QfeLmat(iocn) = tr.QfeL;
    if ~do.qssa
        Jbmat = 0*M3d+NaN; Jbmat(iocn) = Jb;
        Jfmat = 0*M3d+NaN; Jfmat(iocn) = Jf;
    end
end

%% O2
if do.O2
    O2mat = 0*M3d+NaN; O2mat(iocn) = tr.O2;
    Jo2mat = 0*M3d+NaN; Jo2mat(iocn) = Jo2;
    Fo2mat = 0*M3d+NaN; Fo2mat(iocn) = Fo2;
    ResPmat = 0*M3d+NaN; ResPmat(iocn) = res_oxid_p;
    WC_off = 104*(res_oxid_p);
    WC_off_mat = M3d*0+NaN; WC_off_mat(iocn) = WC_off;
    WC_off_prof = squeeze(sum(sum(WC_off_mat.*VOL,1,'omitnan'),2,'omitnan'))*14/1e15;
    WC_off_tot = 14*(WC_off'*dv)/1e15;

    % volume distribution
    bin_width = 2; % in uM O2
    o2bin_bot = [0:bin_width:300];
    o2bin_top = [o2bin_bot(2:end) 300+bin_width];
    o2bin_cen = (o2bin_top + o2bin_top)/2;
    for i = 1:length(o2bin_bot);
        io2 = find(tr.O2>=o2bin_bot(i) & tr.O2<o2bin_top(i));
        vo2(i) = sum(dv(io2));
        io2_obs = find(obs.o2>=o2bin_bot(i) & obs.o2<o2bin_top(i));
        vo2_obs(i) = sum(dv(io2_obs));
    end


end

%% Fe
if do.Fe
    Femat = 0*M3d+NaN; Femat(iocn) = tr.Fe;
    Jfemat = 0*M3d+NaN; Jfemat(iocn) = Jfe;
end

%% Diaz Nutrient limitation
if do.diaz
    Plim = tr.P./(tr.P + KP);
    if do.Fe
        Felim = tr.Fe./(tr.Fe + KFe_diaz);
    else
        Felim = Plim*0 + 1;
    end
    iFe = find(Plim>Felim);
    Lim = ones(m,1);
    Lim(iFe) = 2;
    Lim2 = Lim;
    Lim2(tr.F<1e-10)=0;
    Lim_mat = M3d*0+NaN; Lim2_mat = Lim_mat;
    Lim_mat(iocn) = Lim; Lim2_mat(iocn) = Lim2;
end

%% N:P ratio

if do.var_n2p
    ro_comp = sum(((prod_f(isub).*dv(isub))/(prod_f(isub)'*dv(isub))).*ro(isub));
    vmat = M3d*NaN;
    vmat(iocn) = ro;
end


