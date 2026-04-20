%% Explicit bio


% nutrient limitation
gammaC = tr.P./(tr.P + kpS);
gammaRd = min( gammaC, tr.Fe./(tr.Fe + KFe_diaz) );
gammaC = min( gammaC, tr.Fe./(tr.Fe + kfeS) );
gammaC = min( gammaC, tr.N./(tr.N + knS) );

gammaD = tr.P./(tr.P + kpL);
gammaD = min( gammaD, tr.Fe./(tr.Fe + kfeL) );
gammaD = min( gammaD, tr.N./(tr.N + knL) );
gammaR = (gammaC.*tr.fS)+(gammaD.*tr.fL);

%%
if do.coupled_mub
    mu_b = ((tr.muS.*(tr.fS))+(tr.muL.*(tr.fL)));
    gammaRd = tr.Fe./(tr.Fe + KFe_diaz);
    gammaRd = min( gammaRd, tr.P./(tr.P + KP) );
    gammaC = tr.P./(tr.P + kpS);
    gammaC = min( gammaC, tr.Fe./(tr.Fe + kfeS) );
   prod_b = ((tr.muS.*tr.fS)+(tr.muL.*tr.fL)).*tr.B;%.*gammaR;%.*tlim;
   %mupS = (outS.Vp./tr.QpS).*dpery; mupL = (outL.VpE./tr.QpL).*dpery;
    %munS = (outS.Vn./tr.QnS).*dpery; munL = (outL.Vn./tr.QnL).*dpery;
   %mupS(mupS<0) = 0; mupL(mupL<0) = 0;
   %prod_b = ((mupS.*tr.fS)+(mupL.*tr.fL)).*tr.B;
  % prod_b = (outS.VpE.*tr.fS)+(outL.VpE.*tr.fL);
    prod_f = mu_f.*gammaRd.*tr.F;
else
    prod_f = mu_f.*gammaRd.*tr.F;
    prod_b = mu_b.*gammaR.*tr.B;
end

% compute productivity
prod_b(idif) = 0;
prod_f(idif) = 0;
prod_b(tr.P<=0)=0; prod_b(tr.N<=0)=0; prod_b(tr.Fe<=0)=0;
prod_f(tr.P<=0)=0; prod_f(tr.Fe<=0)=0;



% compute mortality
mort_b = (m1*tr.B + m2*(tr.B+tr.F).*tr.B);
mort_b(tr.B<1e-8)=0;
mort_f = (m1*tr.F + m2*(tr.B+tr.F).*tr.F);
mort_f(tr.F<1e-8)=0;


% SMS
Jb = prod_b - mort_b;
Jf = prod_f - mort_f;


% integrate
% prod_n = production in terms of n

prod_p = prod_b + prod_f;
%prod_n = (outS.VnE.*tr.fS)+(outL.VnE.*tr.fL);
prod_n = prod_b.*ro; % ro: N:P ratio (16:1) - > converting from P to N
S_nfix = prod_f.*rf; % N:P ratio of 50:1

rfe2p2 = (rfe2p_c.*tr.fS)+(rfe2p_d.*tr.fL);
prod_fe = prod_b.*rfe2p2 + prod_f.*rfe2p_diaz;
%prod_fe = prod_b.*((tr.QfeS./tr.QpS.*tr.fs)+(tr.QfeL./tr.Qp.*tr.fs))
mort_p = mort_b + mort_f;
mort_n = mort_b.*ro + mort_f.*rf;
mort_fe = mort_b.*rfe2p2 + mort_f.*rfe2p_diaz;





