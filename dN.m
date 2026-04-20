%% dN

% bio

if ~do.POM;
    exp_n = ef.*mort_n;
    remin_n = Qrem*(exp_n);
    regen_n = (1 - ef - phin).*mort_n;

else % do.POM
    remin_n = tr.pon/taupom;
    regen_n = (1 - phie - phin)*mort_n;
    Jpon = Qsink*tr.pon + phie*mort_n - tr.pon/taupom;
end

Jn = -prod_n + remin_n + regen_n + tr.don/taun;
Jdon = phin*mort_n - tr.don/taun;


% water column denit
if ~ do.O2;
    S_wc = rden*iwc.*remin_p;
else % do.O2
    S_wc = rden*res_oxid_p;
end

% sediment denit
FC = ised.*(106*(phie*mort_p + remin_p).*Zbox*conv);
S_sed = (rden/106)*(exp(-0.9543+0.7662*log(FC)-0.2350*log(FC).*log(FC))/conv)./Zbox;

if do.fixed_denit
    if sum(S_wc)>0 & sum(S_sed)>0;
        % scale to constant value
        cwc = Fwc*Dtot/(14*S_wc'*dv/1e15);
        csed = (1-Fwc)*Dtot/(14*S_sed'*dv/1e15);
        S_wc = cwc*S_wc;
        S_sed = csed*S_sed;
    end
else
    % scale using pre-determined factors
    S_wc = S_wc*wc_scale;
    S_sed = GS_fac.*S_sed*sed_scale;
end

% sum
S_den = do.denit*(S_wc + S_sed);




