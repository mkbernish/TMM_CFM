%% Explicit bio

% compute productivity
gammaR = tr.P./(tr.P + KP); % P limitation
gammaRd = gammaR; %diazotrophs
if do.Fe
    gammaR = min( gammaR, tr.Fe./(tr.Fe + KFe) ); % which one is more
    % limiting, P or Fe
    gammaRd = min( gammaRd, tr.Fe./(tr.Fe + KFe_diaz) ); % same as above for diazos
end
if do.N
    gammaR = min( gammaR, tr.N./(tr.N + KN) );
    % diazos aren't limited by N
end
prod_b = (tr.muS.*tr.fS.*tr.B)+(tr.muL.*tr.fL.*tr.B);
prod_b2 = mu_b.*gammaR.*tr.B; 
prod_b(idif) = 0;
prod_f = mu_f.*gammaRd.*tr.F; 
prod_f(idif) = 0;
prod_b(tr.P<=0)=0; prod_b(tr.N<=0)=0; prod_b(tr.Fe<=0)=0;
prod_f(tr.P<=0)=0; prod_f(tr.Fe<=0)=0;


% compute mortality
mort_b = m1*tr.B + m2*(tr.B+tr.F).*tr.B;
mort_b(tr.B<1e-7)=0;
mort_f = m1*tr.F + m2*(tr.B+tr.F).*tr.F;
mort_f(tr.F<1e-7)=0;


% SMS
Jb = prod_b - mort_b;
Jf = prod_f - mort_f;


% integrate
% prod_n = prod in terms of n
prod_p = prod_b + prod_f;
prod_n = prod_b.*ro; % ro: N:P ratio (16:1) - > converting from P to N
S_nfix = prod_f*rf; % N:P ratio of 50:1
prod_fe = prod_b*rfe2p + prod_f*rfe2p_diaz;

mort_p = mort_b + mort_f;
mort_n = mort_b.*ro + mort_f.*rf;
mort_fe = mort_b*rfe2p + mort_f*rfe2p_diaz;





