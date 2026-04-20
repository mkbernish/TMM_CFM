%% Setup Ninv timestepper %%
%% hovmol plots of various things for tracking / troubleshooting
clear all;
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\"))
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
%close all;
sperd=86400; dpery = 365; spery = sperd*dpery; % time conversions

%% paths

%basepath = '/Users/cdeutsch/Dropbox/projects/TMM/';
basepath = 'C:/Users/mathp/OneDrive/Documents/Papers/paper4/tmm_ss';
mfilepath = "C:/Users/mathp/OneDrive/Documents/Papers/paper4/tmm_ss";

%addpath(genpath(basepath));
addpath(genpath(mfilepath));

mod_mat = '45x90x24';
%mod_mat = '90x180x24';

base_path = [basepath '/Ninv_Fe_v2/'];
model_path = [basepath '/Ninv_Fe_v2/'];
cd(base_path)
%% define run

% control parameters
Dtot = 150;                         % Global denit. rate [TgN/yr]
ro_cons = 16.0;                     % constant n2p ratio
ro_di = 8.0;                        % n2p in non-fixer region
ro_st = 20;                         % n2p for subtropical (picoplankton)
Fe_fac = 1;
zprod = 75;                         % compensation depth
%zprod = 100;                         % compensation depth
% restart
restart = 0;
load_initial = 0;
continue_ts = 0;
initial_file = ['tr_all_DOM_Dfix_150_diazFe_25']; % WARN: CHECK FE_FAC BELOW!!
if restart & load_initial;
    disp(['Cannot do restart AND load_initial']); return; end

% tracers
do.P = 1;
do.O2 = 1;
do.N = 1;
do.Fe = 1;
do.phyto = 1;
do.diaz = 1;
do.DOM = 1; % Dissolved organic matter
do.POM = 1; % Particulate organic matter
if do.diaz & ~do.phyto;
    disp(['Cannot do diaz without phyto']); return; end
if do.phyto & ~do.P;
    disp(['Cannot do phyto without P']); return; end

% processes
do.glacial_O2 = 0;                          % glacial(1 or 5) or modern(0) O2?
do.glacial_Fe = 0;                          % glacial(1) or modern(0) Fe?
do.glacial_SED = 0;                         % glacial(1=weak, 2=strong) or modern(0) sed denit?
do.qssa = 0;                                % QSSA phytoplankton
do.Fe_sms = 1;                              % Fe sources and sinks (closed system if 0)
do.denit = 1;                               % N sources and sinks (closed system if 0)
do.fixed_denit = 0;                         % den
% it feedback if 0
do.var_n2p = 1;                             % variable N:P?
do.spec_rost = 0;                           % specify subtropical N:P?
do.mix_plankton = 0;
do.cons_KN = 1;
do.coupled_fdi = 1;
do.coupled_mub = 1;
do.hov = 1;
do.mov = 0; % 1 if making movie;
do.cons_fe = 0;
do.surf_denit = 0;
do.above_l = 0;

%%
% check compatibility
if do.denit & ~do.diaz;
    disp('WARN: denit but no fixation, N budget will not balance');
end
if do.diaz & ~do.denit;
    disp('WARN: fixation but no denit, N budget will not balance');
end
if do.phyto & do.mix_plankton;
    disp('ERROR Phytoplankton mixing not yet implemented');
    return;
end
%% dfe climatology
load("C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_ss\Ninv_Fe_v2\remapped_dfe.mat")
dfe_remap2 = dfe_remap./1000;
dfe_remap2(38:45,:,:) = 2e-5;
dfe_remap2(isnan(dfe_remap))=nan;

%% get transport matrix & grid

matrix_path = [model_path '/Matrices/grid_matrices'];
Amatrix = 'A_sol2'; % transport matrix to use

load_grid_matrices;
%A = A.*0.1;
%% setup time

ytot = 4;                                % length of sim in year
if m==191169; dt = 1/2000;                  % time step in years
elseif m==45503; dt = 1/700; end            % need short timestep for 2deg model
% shortening timestep
n = (ytot/dt);                              % number of steps
rec = 1;                                   % output frequency in years
if n >10000 && do.mov
    disp('uhhh not happening');
end
if do.mov
    if do.coupled_fdi
        fdi_str = ['cfdi_'];
    else
        fdi_str = ['ofdi_'];
    end

    if do.coupled_mub
        mub_str = ['_cmub'];
    else
        mub_str = ['_omub'];
    end
    mstr = [fdi_str mub_str];
    % mkdir(['movie\'])
    mkdir([ mstr '\N\'])
    mkdir([ mstr '\P\'])
    mkdir([ mstr '\ro\'])
    mkdir([ mstr '\prodn\'])
    mkdir([ mstr '\mu_di\'])
end
%% get data

% WOA 2005
load([model_path, '/Data/WOA05_monthly'], 'obs', '*obs');
%obs.T = obs.T+2;
% Dust (Mahowald 2006)
load(fullfile(model_path, '/Data/Dust_mod_glac'));

% Irradience
par_frac = 0.45; % fraction of solar PAR
k_sw = 0.04; % Attenuation coefficient for PAR (1/depth)
load([model_path, '/Data/isccp_swr'],'*ann','*mon');
Iz = par_frac*repmat(swr_ann,[1 1 length(grid.zt)]).*exp(-k_sw*grid.ZT3d); obs.I = Iz(iocn);

% Oxygen saturation and piston velocity
if do.glacial_O2;
    if do.glacial_O2==1;
        load(fullfile(model_path, '/Data/o2_gasex_glacial'));
    elseif do.glacial_O2==5;
        load(fullfile(model_path, '/Data/o2_gasex_glacial5')); end
else load(fullfile(model_path, '/Data/o2_gasex')); end
obs.kw_o2 = zeros(m,1); obs.kw_o2(ivsurf) = kw_o2(isurf);
obs.o2sat = zeros(m,1); obs.o2sat(ivsurf) = o2sat(isurf);

%% make Q matrices

% make remin matrix (martin curve)
if ~do.POM;
    phie = 0.1;                                         % e-ratio (fraction of NPP to export)
    alpha = 0.858;                                      % shape of martin curve
    [Qneg,Qrem] = qbio(grid,M3d,alpha,kprod);
    Qneg_diag = diag(Qneg);
    Inex = find(Qneg_diag==0);                          % no export (bottom within surface)
    ef = zeros(m,1); ef(isub) = phie; ef(Inex) = 0;     % corrected export ratio
    clear Qneg Qneg_diag

    % make sinking matrix (explicit POM)
else
    phie = 0.1;
    w0 = 3*dpery;                                         % initial sinking rate
    wsink = (1:m)'*0 + w0;
    Qsink = qsink(m,grid,M3d,wsink);
end

save([base_path 'output/' mod_mat '/default_matrices'],'A*','Q*');


%% set up denitrification

% constants and scaling
rden = 350;
Fwc = 1/3;
Zbox = grid.DZT3d(iocn);
conv = 1000/(1e4*365.25);
wc_scale = 0.25;
sed_scale = 0.5;

% find bottom boxes
IBOT = repmat(max(i3d,[],3),[1 1 length(grid.zt)]);
ised = zeros(m,1); ised(IBOT(iocn)==i3d(iocn)) = 1;

% find low observed o2
iwc = zeros(m,1); iwc(obs.o2<10) = 1;

% glacial sediment denit
GS_cont = 0.25;
GS_fac = ones(m,1);
if do.glacial_SED==1;
    GS_fac(icont) = 0.25; end
if do.glacial_SED==2;
    GS_fac = GS_fac*0.4; end

%% setup dust deposition

%fFe = 0.035;
fFe = 0.01;

% modern
if ~ do.glacial_Fe
    dust_dep = zeros(m,1); dust_dep(ivsurf) = dust(isurf);
    fe_dep = fFe*dust_dep;
    % glacial
else
    dust_dep = zeros(m,1); dust_dep(ivsurf) = dust_glac(isurf);
    fe_dep = fFe*dust_dep;
end


%% bio parameters

% growth and mort
mumax = 1*dpery;                  % max growth rate of non-fixer
mumax2 = 1*dpery;
fix_cost = 0.85;                    % 'cost' of fixation
KP = 0.25;                           % P half saturation [uM]
m1 = 0.1*dpery;                     % linear mortality
m2 = 10*dpery;                    % quadratic mortality
T_nofix = 10.0;                      % threshold temperature for N2 fixation

% ratios
rf = 30;                            % fixer N/P ratio
KSi = 6;                          % phyto half-sat for Silicate

ro2p = -150;
% r = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_cfdi__cmub_varFe_5_reg16.mat');

rfe2p = 2e-4;
rfe2p_c = 1.8e-4;
rfe2p_d = 1.2e-4;
KN_cons = ro_cons*KP;

% limitation
To = 30;                            % max temperature
k_C = 0;                          % temperature sensitivity
k_D = 0;                          % temperature sensitivity
k_diaz = 0.025;
KI = 20;                             % I half saturation
Kdust = 0.25;                       % 'half-deposition' for fixers

% dom and pom
phin = do.DOM*0.1;                  % fraction converted to DON
phip = do.DOM*0.2;                  % fraction converted to DOP
% phip has to be higher than phin
taun = 0.8;                           % lifetime of DON in year
taup = 0.5;                           % lifetime of DOP in year
% taup has to be smaller than taun -
%taupom = 1/(0.033*dpery);           % lifetime of POM in year
taupom = 1/(0.02*dpery);           % lifetime of POM in year
% Fe
alpha_fe = 0.01;                   % Solubility
Ligand = 1/1000;                    % Fe binding ligand concentration
beta_fe = 2e5;                      % Fe ligand binding concentration (uM^-1)
ksc = 1.1e-3*dpery;                 % Fe scavenging (year^-1)

% Misc.
bio_tau = 1/12;                     % nutrient damping timescale
eps_lim = 1e-6;                     % minimum tracer conc.
o2_crit = 7;                        % critical oxygen
KO2 = 0;                            % remineralization O2 limitation (half-sat)

%% prescribe model name

make_model_name;
nmod = [nmod '_' num2str(ytot) '_reg19']

%% bio limitation

ilim = (1-exp(-(obs.I/KI)));
tlimC = exp(k_C*(obs.T-To)); tlimC(obs.T>To) = 1;
tlimD = exp(k_D*(obs.T-To)); tlimD(obs.T>To) = 1;
tlim_diaz = exp(k_diaz*(obs.T-To)); tlim_diaz(obs.T>To) = 1;
mu_b_const =  mumax2.*ilim.*tlimD;
mu_b = mumax.*ilim.*tlimD;  % growth rate, non-fixing plankton
mu_f = mumax.*fix_cost.*ilim.*tlim_diaz; % growth rate, N-fixing plankton
mu_f(obs.T < T_nofix) = 0;

%% N/P ratio

if ~do.var_n2p
    ro = ro_cons;

else % variable N:P - tied to the spatial dist. of observed silicate
    fdi = obs.Si./(obs.Si+KSi); % fraction of community diatoms
    ro_ndi = 16+(16-ro_di); %24
    ro = fdi*ro_di + (1-fdi)*ro_ndi;
end

if ~do.cons_KN
    KN = ro*KP;
else
    KN = KN_cons;
end

if do.cons_fe
    load('remapped_dfe.mat')
    fev = dfe_remap(iocn);
end

%% initial conditions
v_fe = dfe_remap2(iocn);
if restart
    try load(['output/' mod_mat '/'  nmod],'tr','ro*');
        if ~do.cons_KN; KN = ro*KP; end
        disp([nmod ' loaded initial conditions from restart file']);
        if continue_ts;
            load(['output/' mod_mat '/' nmod],'*inv','source*','sink*','time');
            TS_end = length(P_inv);
            year0 = time(end);
            disp('Loaded timeseries'); end
    catch
        disp(['Error: couldnt find restart file']); return
    end


elseif load_initial
    try load(['output/' mod_mat '/' initial_file],'tr','ro');
        if ~do.cons_KN; KN = ro*KP; end
        disp([nmod ' loaded initial conditions from file ' initial_file]);
        year0 = 0;
        if continue_ts;
            load(['output/' mod_mat '/' initial_file],'*inv','source*','sink*','time');
            TS_end = length(P_inv);
            year0 = time(end);
            disp('Loaded timeseries'); end
    catch
        disp(['Error: couldnt find initail conditions file']); return
    end

else % use observed or uniform initial conditions
    % tr - anything that gets transported or transformed
    tr.P = obs.P;
    tr.N = obs.N;
    tr.O2 = obs.o2;
    % if do.cons_fe
    %     tr.Fe = fev./1000;
    % else
    %     tr.Fe = fe_dep;
    % end
    % tr.Fe = zeros(m,1) + 6e-4;
    tr.Fe = v_fe;
    tr.dop = zeros(m,1);
    tr.don = zeros(m,1);
    tr.pop = zeros(m,1);
    tr.pon = zeros(m,1);
    tr.pofe = zeros(m,1);
    tr.B = zeros(m,1) + 0.01;
    if do.phyto; tr.B(isub)=1e-2.*tlim_diaz(isub); end
    % B is in mmol P m-3 (Phytoplankton biomass)
    tr.QnS = zeros(m,1)+(0.23); % in mol N per mol C
    tr.QpS = zeros(m,1)+(0.02); % in mol N per mol C
    tr.QfeS = zeros(m,1)+(0.00001); % in mol Fe per mol C
    tr.QnL = zeros(m,1)+(0.23); % in mol N per mol C
    tr.QpL = zeros(m,1)+(0.024); % in mol N per mol C
    tr.QfeL = zeros(m,1)+(0.00001); % in mol Fe per mol C
    tr.muS= zeros(m,1)+(0.4.*dpery);
    tr.muL = zeros(m,1)+(0.5.*dpery);
    tr.fL= obs.Si./(obs.Si+KSi);
    tr.fS = 1-tr.fL;
    % tr.ro = zeros(m,1)+6.625;
    tr.F = zeros(m,1); if do.diaz; tr.F(isub)=1e-2.*tlim_diaz(isub); end
    disp([nmod ' used observed or uniform initial conditions']);
    year0 = 0;
    if ~do.cons_KN
        KN = ro*KP;
    end
end

if ~exist('TS_end');
    TS_end = 0;
    year0 = 0; end

disp(['Simulation starting at year ' num2str(year0)]);

% timestepper initial conditions
% g values keep track of time rate of change
efac = 0.1;
gPm1 = zeros(m,1); % tendency at time step i-1
gNm1 = zeros(m,1);
gFem1 = zeros(m,1);
gOm1 = zeros(m,1);
gBm1 = zeros(m,1);
gFm1 = zeros(m,1);
gDPm1 = zeros(m,1);
gDNm1 = zeros(m,1);
gPPm1 = zeros(m,1);
gPNm1 = zeros(m,1);
gPFem1 = zeros(m,1);
gQnm1S = zeros(m,1); % small
gQpm1S = zeros(m,1);
gQfem1S = zeros(m,1);
gQnm1L = zeros(m,1); % large
gQpm1L = zeros(m,1);
gQfem1L = zeros(m,1);
outS.Vp = zeros(m,1);
outL.Vp = zeros(m,1);
outS.VpE = zeros(m,1);
outL.VpE = zeros(m,1);
bigtso = nan([n,82]);
bigtna = nan([n,82]);
bigtpg = nan([n,82]);


%% 2022
% PmaxS = 0.0035;
% %PmaxS = 0.1;
% knS=0.8; kpS=0.2;
% kfeS=0.1e-4;
% vmaxnS = 0.1; vmaxpS = 0.012; vmaxfeS = 0.15e-4;
% QpmaxS = 1.4e-2;A_Fe_phoS = 0.28e-2; EaS = 67800;
%
% PmaxL = 0.006;
% %PmaxL = 0.2;
% knL =1.4; kpL =0.38;
% kfeL =0.2e-4;
% vmaxnL =0.2; vmaxpL =0.037; vmaxfeL =0.3e-4;
% QpmaxL = 4e-2;A_Fe_phoL = 0.12e-2; EaL = 33700;
%
% rfe2p_diaz = rfe2p_c*1.01;
% KFe_diaz = 4e-4;

%% lowering everything
%% 2022
PmaxS = 0.003;
knS=0.7; kpS=0.35;
kfeS=3e-4;
vmaxnS = 0.4; vmaxpS = 0.015; vmaxfeS =0.6e-4;
QpmaxS = 1.8e-2;A_Fe_phoS = 6e-3; EaS = 70000;

PmaxL = 0.0039;
knL =0.9; kpL =0.55;
kfeL =4e-4;
vmaxnL =0.5; vmaxpL =0.022; vmaxfeL =1.2e-4;
QpmaxL =2.4e-2;A_Fe_phoL = 5e-3; EaL = 30000;

rfe2p_diaz = rfe2p_c*1.02;
KFe_diaz = 3.4e-4;
KP = 0.35;

% PmaxS = 0.0032;
% knS=0.4; kpS=0.5;
% kfeS=3e-4;
% vmaxnS = 0.4; vmaxpS = 0.03; vmaxfeS =1e-4;
% QpmaxS = 2e-2;A_Fe_phoS = 5.8e-3; EaS = 80000;
% 
% PmaxL = 0.0039;
% knL =0.55; kpL =0.65;
% kfeL =3.8e-4;
% vmaxnL =0.55; vmaxpL =0.04; vmaxfeL =1.6e-4;
% QpmaxL =2.8e-2;A_Fe_phoL = 4.8e-3; EaL = 30000;
% 
% rfe2p_diaz = rfe2p_c*1.02;
% KFe_diaz = 3e-4;
% KP = 0.48;

%% differentiating
% PmaxS = 0.0032;
% knS=0.4; kpS=0.1;
% kfeS=2e-4;
% vmaxnS = 0.4; vmaxpS = 0.03; vmaxfeS = 0.8e-4;
% QpmaxS = 2.4e-2;A_Fe_phoS = 5.8e-3; EaS = 70000;
%
% PmaxL = 0.0039;
% knL =0.6; kpL =0.2;
% kfeL =3.5e-4;
% vmaxnL =0.5; vmaxpL =0.05; vmaxfeL =1.4e-4;
% QpmaxL = 3e-2;A_Fe_phoL = 4.8e-3; EaL = 33700;
%
% rfe2p_diaz = rfe2p_c*1.02;
% KFe_diaz = 2e-4;
% KP = kpS;

%% loading old params
% PmaxS = r.PmaxS;
% knS=r.knS; kpS=r.kpS;
% kfeS=r.kfeS;
% vmaxnS = r.vmaxnS; vmaxpS = r.vmaxpS; vmaxfeS = r.vmaxfeS;
% QpmaxS = r.QpmaxS;A_Fe_phoS = 5.83e-3; EaS = r.EaS;
%
% PmaxL = r.PmaxL;
% knL =r.knL; kpL =r.kpL;
% kfeL =r.kfeL;
% vmaxnL =r.vmaxnL; vmaxpL =r.vmaxpL; vmaxfeL =r.vmaxfeL;
% QpmaxL = r.QpmaxL;A_Fe_phoL = 4.8e-3; EaL = r.EaL;
%
% rfe2p_diaz = rfe2p_c*1.02;
% KFe_diaz = 1.08e-4;


%% send to timestepper

tic
time_step
t = toc;

disp([nmod ' finished ' num2str(ytot) ' year simulation in ' num2str(t/60) ' minutes']);

bigtso2 = array2table(bigtso,'VariableNames',{'N','P','Fe','mundi','mudi','B',...
    'prodB','prodN','prodP','regenP','reminP','regenN','reminN','reminFe','regenFe','DON','DOP',...
    'mortb','mortp','mortn','mub','ro',...
    'QnS','QnL','QpS','QpL','QfeS','QfeL',...
    'Jn','gN','prod_fe','F','mu_f',...
    'muNS','muPS','muFeS','muCS','mumaxS',...
    'muNL','muPL','muFeL','muCL','mumaxL',...
    'fracPhoS','fracBioS','fracStoS','fracStrS','mort_f',...
    'S_den','prod_f','S_nfix',...
    'NbioS','NphoS','NstoS','NchlS','NdnaS','NessS','NprotS','NrnaS','Nprot_bS','Nprot_pS','Nprot_cS','VnS','Npho2S',...
    'NbioL','NphoL','NstoL','NchlL','NdnaL','NessL','NprotL','NrnaL','Nprot_bL','Nprot_pL','Nprot_cL','VnL','Npho2L',...
    'aN','VnES','VnEL','VpS','VpL'});
bigtna2 = array2table(bigtna,'VariableNames',{'N','P','Fe','mundi','mudi','B',...
    'prodB','prodN','prodP','regenP','reminP','regenN','reminN','reminFe','regenFe','DON','DOP',...
    'mortb','mortp','mortn','mub','ro',...
    'QnS','QnL','QpS','QpL','QfeS','QfeL',...
    'Jn','gN','prod_fe','F','mu_f',...
    'muNS','muPS','muFeS','muCS','mumaxS',...
    'muNL','muPL','muFeL','muCL','mumaxL',...
    'fracPhoS','fracBioS','fracStoS','fracStrS','mort_f',...
    'S_den','prod_f','S_nfix',...
    'NbioS','NphoS','NstoS','NchlS','NdnaS','NessS','NprotS','NrnaS','Nprot_bS','Nprot_pS','Nprot_cS','VnS','Npho2S',...
    'NbioL','NphoL','NstoL','NchlL','NdnaL','NessL','NprotL','NrnaL','Nprot_bL','Nprot_pL','Nprot_cL','VnL','Npho2L',...
    'aN','VnES','VnEL','VpS','VpL'});

bigtpg2 = array2table(bigtpg,'VariableNames',{'N','P','Fe','mundi','mudi','B',...
    'prodB','prodN','prodP','regenP','reminP','regenN','reminN','reminFe','regenFe','DON','DOP',...
    'mortb','mortp','mortn','mub','ro',...
    'QnS','QnL','QpS','QpL','QfeS','QfeL',...
    'Jn','gN','prod_fe','F','mu_f',...
    'muNS','muPS','muFeS','muCS','mumaxS',...
    'muNL','muPL','muFeL','muCL','mumaxL',...
    'fracPhoS','fracBioS','fracStoS','fracStrS','mort_f',...
    'S_den','prod_f','S_nfix',...
    'NbioS','NphoS','NstoS','NchlS','NdnaS','NessS','NprotS','NrnaS','Nprot_bS','Nprot_pS','Nprot_cS','VnS','Npho2S',...
    'NbioL','NphoL','NstoL','NchlL','NdnaL','NessL','NprotL','NrnaL','Nprot_bL','Nprot_pL','Nprot_cL','VnL','Npho2L',...
    'aN','VnES','VnEL','VpS','VpL'});


%
% outS.Nbio(ix2),outS.Npho(ix2),outS.Nstor(ix2),outS.Nchl(ix2),outS.Ndna(ix2),outS.Ness(ix2),outS.Nprotein(ix2),...
% outS.Nrna(ix2),outS.Nprot_biosynth(ix2),outS.Nprot_photo(ix2),outS.Nconst_prot(ix2),outS.Vn(ix2)...
%          outL.Nbio(ix2),outL.Npho(ix2),outL.Nstor(ix2),outL.Nchl(ix2),outL.Ndna(ix2),outL.Ness(ix2),outL.Nprotein(ix2),...
% outL.Nrna(ix2),outL.Nprot_biosynth(ix2),outL.Nprot_photo(ix2),outL.Nconst_prot(ix2),outL.Vn(ix2)

%%

calc_results
clear A*
save([base_path 'output/' mod_mat '/' nmod],'*');
%
%%
plt.tst = 0;
if plt.tst
    fs = 9;
    ww = 0*M3d+NaN; ww(iocn) =ro;
    cmap = getPyPlot_cMap('jet', 24,'pyCmd','py');
    figure('Position',[100 100 900 850])
    [ha,pos]=tight_subplot(4,2,[0.02 0.02],[.04 .08],[.04 .05]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),linspace(6,30,24),'linecolor','none')
    shading interp;
    clim([6 30])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);

    axes(ha(2))
    ww = 0*M3d+NaN; ww(iocn) =ro;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)),20,'linecolor','none')
    shading interp;
    clim([6 30])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);

    axes(ha(3))
    ww = 0*M3d+NaN; ww(iocn) =tr.QnS;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),12,'linecolor','none')
    shading interp;
    clim([0 0.24])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);

    axes(ha(4))
    ww = 0*M3d+NaN; ww(iocn) =outS.Qn_max-tr.QnS;
    Nstor_max = 2e-15./outS.Qc;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),12,'linecolor','none')
    shading interp;
   % clim([0 0.24])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);
    clim([0 0.024])
    

    axes(ha(5))
    ww = 0*M3d+NaN; ww(iocn) =tr.QpS;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),12,'linecolor','none')
    shading interp;
    clim([0 0.024])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);

    axes(ha(6))
    ww = 0*M3d+NaN; ww(iocn) =tr.QpS;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)),12,'linecolor','none')
    shading interp;
    clim([0 0.024])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);

    axes(ha(7))
    ww = 0*M3d+NaN; ww(iocn) =tr.ro_ndiv;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),24,'linecolor','none')
    shading interp;
    clim([6 30])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);

    axes(ha(8))
    ww = 0*M3d+NaN; ww(iocn) =tr.ro_ndiv;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)),24,'linecolor','none')
    shading interp;
    clim([6 30])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);
end
%% lat NP
plt.uh = 1;
if plt.uh
    gr = grid;
    tt1 = readtable("C:\Users\mathp\OneDrive\Documents\Papers\paper4\real_data\poc_pondata_2025.xlsx");
    bigt = tt1;
    cmap = flipud(getPyPlot_cMap('RdYlBu', 256,'pyCmd','py'));
    fs = 11;
    bigtsurf = bigt(bigt.Depth<=75,:);
    bigt2 = bigtsurf;
    bigt2(bigtsurf.Longitude>0 & bigtsurf.Longitude<60&...
        bigtsurf.Latitude>30 &bigtsurf.Latitude<60,:) = [];
    [u,~,j]=unique(bigt2.Latitude,'rows','stable');
    latpon=[  accumarray(j,bigt2.C_N,[],@mean)  , u];
    latpop=[  accumarray(j,bigt2.C_P,[],@mean)  , u];
    latnp=[  accumarray(j,bigt2.N_P,[],@mean)  , u];
    figure('Position',[950 200 750 550])
    [ha,pos]=tight_subplot(1,1,[0.03 0.03],[.24 .2],[.1 .04]);

    nps = ro;

    romatc = 0*M3d+NaN; romatc(iocn) = nps;
    err = std(mean(romatc(:,:,1:2),3,'omitnan'),[],2,'omitnan');
    nps_atl = romatc;nps_atl(MSKS.ATL==0)=nan;
    nps_pac = romatc;nps_pac(MSKS.PAC==0)=nan;
    nps_ind = romatc;nps_ind(MSKS.IND==0)=nan;

    axes(ha(1))
    scatter(latpon(:,2),latnp(:,1),30 ...
        ,'k','filled','MarkerFaceAlpha',0.1,'markeredgecolor','k','MarkerEdgeAlpha',0.3)
    hold on;
    plot(gr.yt,mean(nps_atl(:,:,1:2),[2 3],'omitnan'),'b-','linewidth',1.2)
    hold on;
    plot(gr.yt,mean(nps_pac(:,:,1:2),[2 3],'omitnan'),'r-','linewidth',1.2)
    plot(gr.yt,mean(nps_ind(:,:,1:2),[2 3],'omitnan'),'g-','linewidth',1.2)
    errorbar(gr.yt,mean(romatc(:,:,1:2),[2 3],'omitnan'),err,'m-','linewidth',0.8);

    l1=legend({'Observed (Liu et al., 2025)','ATL','PAC','IND'},'location','northeastoutside',...
        'numcolumns',4);
    l1.Position = [0.28    0.81    0.2771    0.04];
    xlim([-75 75])
    grid on; box on;
    ylim([5 45])
    set(gca,'fontsize',12,'xticklabel',{''})
    ylabel('Planktonic N:P (mol mol^{-1})')
    % set(gca,'yscale','log')
    xticklabels({['60' char(176) 'S'],...
        ['40' char(176) 'S'], ['20' char(176) 'S'],...
        ['0' char(176)], ['20' char(176) 'N'],...
        ['40' char(176) 'N'], ['60' char(176) 'N']})
end

%% one
plt.one = 0;
if plt.one
    fs = 12;
    ww = 0*M3d+NaN; ww(iocn) =ro;
    cmap = flipud(getPyPlot_cMap('RdYlBu', 22,'pyCmd','py'));
    figure('Position',[900 200 600 750])
    [ha,pos]=tight_subplot(2,1,[0.02 0.02],[.04 .08],[.04 .05]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)),linspace(8,30,22),'linecolor','none')
    shading interp;
    colormap(cmap(:,1:3))
    clim([8 30])
    %cmocean('balance',24,'pivot',16)
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120]);
    cb1=colorbar('horizontal','fontsize',fs);

    axes(ha(2))
    ww = 0*M3d+NaN; ww(iocn) =tr.P;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)),50,'linecolor','none')
    shading interp;
    clim([0 1])
    %cmocean('balance',50,'pivot',0.5)
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120]);
    cb2=colorbar('horizontal','fontsize',fs);
end

%% np1
plt.np =1;
if plt.np
    fs = 12;
    ww = 0*M3d+NaN; ww(iocn) =ro;
    cmap = flipud(getPyPlot_cMap('RdYlBu', 256,'pyCmd','py'));
    figure('Position',[100 100 850 700])
    [ha,pos]=tight_subplot(3,2,[0.04 0.04],[.1 .1],[.06 .03]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(8,30,11);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)))
    hold on;shading interp;
    clim([8 32])
   % clim([8 36])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120]);
    cb1=colorbar('fontsize',fs);
    % cb1.Position = [ 0.12    0.77   0.32   0.024];
    cb1.Label.String = 'R_O (mol N (mol P)^{-1})';
    % cb1.Ticks = [8, 12, 16, 20, 24, 28, 32

    axes(ha(2))
    ww2 = 0*M3d+NaN; ww2(iocn) =tr.pon./tr.pop;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww2(:,:,1:2),3,'omitnan'),45,2)))
    hold on;shading interp;
  clim([8 32])
    m_contour(grid.xt-180,grid.yt,...
        (circshift(mean(ww2(:,:,1),3,'omitnan'),45,2)),[0 0],'linecolor','g')
    m_contour(grid.xt-180,grid.yt,...
        (circshift(mean(ww2(:,:,2),3,'omitnan'),45,2)),[0 0],'linecolor','r')
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120]);
    cb2=colorbar('fontsize',fs);
    % cb2.Position = [ 0.59    0.77   0.32   0.024];
    cb2.Label.String = 'Diatom fraction';
    %cb2.Ticks = [0, 0.25,0.5,0.75,1];


    axes(ha(3))
    %  nsto_max = 2.9e-15./outL.Qc;
    % qn_max = (outL.Qn_nostor./outL.Qc)+ nsto_max;
    ww = 0*M3d+NaN; ww(iocn) =tr.N;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)))
    hold on;shading interp;
    m_contour(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),[0 0],'linecolor','g')
    m_contour(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)),[0 0],'linecolor','r')
    clim([0 30])
    %clim([-4 4])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120]);
    cb3=colorbar('fontsize',fs);
    % cb2.Position = [ 0.59    0.77   0.32   0.024];
    cb3.Label.String = 'Diatom fraction';
    %cb2.Ticks = [0, 0.25,0.5,0.75,1];

    axes(ha(4))
    ww = 0*M3d+NaN; ww(iocn) =outS.limType;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)))
    hold on;shading interp;
    clim([1 4])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120]);
    cb3=colorbar('fontsize',fs);
    cb3.Label.String = 'Diatom fraction';


    axes(ha(5))
    ww = 0*M3d+NaN; ww(iocn) =outS.limType;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)))
    hold on;shading interp;
    clim([1 4])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120]);
    cb3=colorbar('fontsize',fs);
    cb3.Label.String = 'Diatom fraction';

    axes(ha(6))
    ww = 0*M3d+NaN; ww(iocn) =tr.B;
  %  ww(ww<1e-9) = nan;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)))
    hold on;shading interp;
    m_contour(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),[0 0],'linecolor','g')
    m_contour(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)),[0 0],'linecolor','r')
    %clim([0 2])
    clim([0 0.02])
   % clim([0 0.04])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120]);
    cb3=colorbar('fontsize',fs);
    cb3.Label.String = 'F';
    dd = 75;


end

%% plt.reg
plt.reg = 0;
if plt.reg
    %   reg = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_fc75_cfdi__cmub_varFe_1_reg19.mat');
    fs = 10;
    z = 1:2;
    ww = 0*M3d+NaN; ww(iocn) =ro;
    cmap = flipud(getPyPlot_cMap('RdYlBu', 256,'pyCmd','py'));
    figure('Position',[100 200 800 600])
    [ha,pos]=tight_subplot(3,2,[0.02 0.02],[.04 .04],[.04 .05]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[1 360]);
    m_contourf(grid.xt,grid.yt,...
        (mean(ww(:,:,z),3,'omitnan')),100,'linecolor','none')
    shading interp;
    colormap(cmap(:,1:3))
    clim([8 30])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    Ml = biomass_npD(end);Ms = biomass_npC(end);
    np_size = Ml.*tr.fL+(Ms.*tr.fS);

    axes(ha(2))
    ww = 0*M3d+NaN; ww(iocn) =np_size;
    m_proj('miller','lat',[-80 80],...
        'lon',[1 360]);
    m_contourf(grid.xt,grid.yt,...
        (mean(ww(:,:,z),3,'omitnan')),100,'linecolor','none')
    shading interp;
    colormap(cmap(30:200,1:3))
    clim([10 24])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);

    np_sig = ro-np_size;
    axes(ha(3))
    ww = 0*M3d+NaN; ww(iocn) =np_sig;
    m_proj('miller','lat',[-80 80],...
        'lon',[1 360]);
    m_contourf(grid.xt,grid.yt,...
        (mean(ww(:,:,z),3,'omitnan')),100,'linecolor','none')
    hold on;
    m_contour(grid.xt,grid.yt,...
        (mean(ww(:,:,z),3,'omitnan')),[0 0],'linecolor','w')
    shading interp;
    colormap(cmap(:,1:3))
    clim([-10 10])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);

    pstor = (outS.Pstor./outS.Qc.*tr.fS)+...
        (outL.Pstor./outL.Qc.*tr.fL);
    pstorn = pstor./(tr.QnS.*tr.fS+(tr.QnL.*tr.fL));

    axes(ha(4))
    ww = 0*M3d+NaN; ww(iocn) =pstorn;
    m_proj('miller','lat',[-80 80],...
        'lon',[1 360]);
    m_pcolor(grid.xt,grid.yt,...
        (mean(ww(:,:,z),3,'omitnan')))
    hold on;
    shading interp;
    colormap(cmap(:,1:3))
    clim([0 0.12])
    % clim([-8 8])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);

    axes(ha(5))
    ww = 0*M3d+NaN; ww(iocn) =1:m;
    m_proj('miller','lat',[-80 80],...
        'lon',[1 360]);
    m_contourf(grid.xt,grid.yt,...
        (mean(ww(:,:,1),3,'omitnan')),100,'linecolor','none')
    hold on;
    shading flat;
    % clim([-10 10])
    colormap(cmap(:,1:3))
    %
    %clim([0 0.5])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);

    axes(ha(6))
    ww = 0*M3d+NaN; ww(iocn) =outS.limType;
    m_proj('miller','lat',[-80 80],...
        'lon',[1 360]);
    m_contourf(grid.xt,grid.yt,...
        (mean(ww(:,:,1),3,'omitnan')),[1 2 3 4],'linecolor','none')
    hold on;
    shading interp;
    colormap(cmap(:,1:3))
    % clim([0 2])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
end

%% Pacific
plt.np2 =0;
if plt.np2
    fs = 9;
    ww = 0*M3d+NaN; ww(iocn) =ro;
    cmap = getPyPlot_cMap('jet', 256,'pyCmd','py');
    figure('Position',[50 100 850 760])
    [ha,pos]=tight_subplot(6,3,[0.02 0.02],[.02 .02],[.02 .02]);
    axes(ha(1))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),17,'linecolor','none')
    hold on;shading interp;
    clim([8 25])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'R_O (z=1)','fontsize',fs)
    %

    axes(ha(2))
    ww = 0*M3d+NaN; ww(iocn) =ro;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),240,2)),17,'linecolor','none')
    hold on;shading interp;
    % clim([8 32])
    clim([8 25])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'R_O (z=2)','fontsize',fs)

    axes(ha(3))
    ww = 0*M3d+NaN; ww(iocn) =ro;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),240,2)),17,'linecolor','none')
    hold on;shading interp;
    clim([8 25])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'R_O (avg)','fontsize',fs)
    %
    axes(ha(4))
    ww = 0*M3d+NaN; ww(iocn) =tr.fL;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),10,'linecolor','none')
    hold on;shading interp;
    clim([0 1])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'f_D (z=1)','fontsize',fs)

    axes(ha(5))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),240,2)),10,'linecolor','none')
    hold on;shading interp;
    clim([0 1])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'f_D (z=2)','fontsize',fs)

    axes(ha(6))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),240,2)),10,'linecolor','none')
    hold on;shading interp;
    clim([0 1])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'f_D (avg)','fontsize',fs)

    axes(ha(7))
    ww = 0*M3d+NaN; ww(iocn) =tr.QnS;
    ww2 = 0*M3d+NaN; ww2(iocn) =outS.Qn_max;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_pcolor(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)))
    hold on;shading interp;
        m_contour(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)))
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'Q_N (z=1)','fontsize',fs)
    clim([0.1 0.24])

    axes(ha(8))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_pcolor(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),240,2)))
    hold on;shading interp;
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'Q_N (z=2)','fontsize',fs)
    clim([0.1 0.24])

    axes(ha(9))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_pcolor(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),240,2)))
    hold on;shading interp;
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    m_text(-220,44,'Q_N (avg)','fontsize',fs)
    clim([0.1 0.24])
    %

    axes(ha(10))
    ww = 0*M3d+NaN; ww(iocn) =(tr.QpS.*tr.fS)+(tr.QpL.*tr.fL);
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'Q_P (z=1)','fontsize',fs)
    clim([0.002 0.012])

    axes(ha(11))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,2),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'Q_P (z=2)','fontsize',fs)
    clim([0.002 0.012])

    axes(ha(12))
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1:2),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'Q_P (avg)','fontsize',fs)
    clim([0.005 0.016])

    axes(ha(13))
    ww = 0*M3d+NaN; ww(iocn) =tr.N;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
      m_contour(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),[1 1],'linecolor','m')
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'N (z=1)','fontsize',fs)
    clim([0 10])

    axes(ha(14))
     ww = 0*M3d+NaN; ww(iocn) =tr.N-obs.N;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),120,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'mod-obs N','fontsize',fs)
    clim([-2 2])

    axes(ha(15))
    ww = 0*M3d+NaN; ww(iocn) =tr.P;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'P (z=1)','fontsize',fs)
    clim([0 1])


    axes(ha(16))
    ww = 0*M3d+NaN; ww(iocn) =tr.ro_ndiv;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'NP_{cyano} (z=1)','fontsize',fs)
    clim([8 28])

    axes(ha(17))
    ww = 0*M3d+NaN; ww(iocn) =tr.ro_div;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'NP_{diat} (z=2)','fontsize',fs)
    clim([8 28])

    axes(ha(18))
    ww = 0*M3d+NaN; ww(iocn) =tr.Fe.*1000;
    m_proj('miller',['lat'],[-40 40],...
        'lon',[-240 -60]);
    m_contourf(grid.xt-240,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),240,2)),20,'linecolor','none')
    hold on;shading interp;
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb3=colorbar('fontsize',fs);
    m_text(-220,44,'Fe (z=1)','fontsize',fs)
    clim([0 1])


end

%% allocation map
plt.alloc = 0;
if plt.alloc
    alloc_type = max([outS.fracPhoto, outS.fracBiosynth, outS.fracStor,...
        outS.fracStruct_const],[],2) ;
    att = nan(size(outS.fracPhoto));
    for ii = 1:length(att)
        if alloc_type(ii) == outS.fracPhoto(ii)
            att(ii) = 1 ;    %photosyn
        elseif alloc_type(ii) == outS.fracBiosynth(ii)
            att(ii) = 3 ;    % bio
        elseif alloc_type(ii) == outS.fracStor(ii)
            att(ii) = 4 ;    %storage
        elseif alloc_type(ii) == outS.fracStruct_const(ii)
            att(ii) = 2 ;    % structual/essential
        end
    end
    cmap = flipud(getPyPlot_cMap('Pastel1', 256,'pyCmd','py'));
    Lim_cmap = [cmap(190,1:3); ... % green for photosynthesis
        cmap(80,1:3); ... % brownish
        cmap(150,1:3); ... % purple
        cmap(240,1:3)] ;%reddish pink
    fs = 10;
    ww = 0*M3d+NaN; ww(iocn) =att;
    figure('Position',[100 100 700 750])
    [ha,pos]=tight_subplot(2,1,[0.04 0.04],[.08 .06],[.04 .12]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(8,30,11);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)))
    hold on;shading flat;
    set(gca, 'CLim', [0.5 4.5],'Colormap', Lim_cmap)
    cbh = colorbar('Ticks', [1, 2, 3, 4], 'TickLabels',...
        {'Photosynthesis', 'Essential/Structural', 'Biosynthesis', 'Storage'});
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120]);


    alloc_type = max([outL.fracPhoto, outL.fracBiosynth, outL.fracStor,...
        outL.fracStruct_const],[],2) ;
    att = nan(size(outL.fracPhoto));
    for ii = 1:length(att)
        if alloc_type(ii) == outL.fracPhoto(ii)
            att(ii) = 1 ;    %photosyn
        elseif alloc_type(ii) == outL.fracBiosynth(ii)
            att(ii) = 3 ;    % bio
        elseif alloc_type(ii) == outL.fracStor(ii)
            att(ii) = 4 ;    %storage
        elseif alloc_type(ii) == outL.fracStruct_const(ii)
            att(ii) = 2 ;    % structual/essential
        end
    end
    cmap = flipud(getPyPlot_cMap('Pastel1', 256,'pyCmd','py'));
    Lim_cmap = [cmap(190,1:3); ... % green for photosynthesis
        cmap(80,1:3); ... % brownish
        cmap(150,1:3); ... % purple
        cmap(240,1:3)] ;%reddish pink
    fs = 10;
    ww = 0*M3d+NaN; ww(iocn) =att;
    axes(ha(2))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(8,30,11);
    m_pcolor(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)))
    hold on;shading flat;
    set(gca, 'CLim', [0.5 4.5],'Colormap', Lim_cmap)
    cbh = colorbar('Ticks', [1, 2, 3, 4], 'TickLabels',...
        {'Photosynthesis', 'Essential/Structural', 'Biosynthesis', 'Storage'});
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120]);
end


%% comparing two
plt.lp=0;
if plt.lp
    xl = [1 3200];
    figure('Position',[50 60 730 950])
    [ha,pos]=tight_subplot(11,2,[0.02 0.04],[.02 .02],[.07 .02]);
   
    axes(ha(1))
    plot(1:n,bigtso2.N,'k-') % bigtso - outside triangle, low NP
    hold on;
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('[N]')
    xlim([xl(1) xl(2)])
    title('Low N:P')
    ylim([0 0.2])

    axes(ha(2))
    plot(1:n,bigtna2.N,'k-')
    hold on;
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
    title('High N:P')
   ylim([0 0.2])


    axes(ha(3))
    plot(1:n,bigtso2.P,'k-')
    hold on;
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('[P]')
    xlim([xl(1) xl(2)])
 ylim([0 1])

    axes(ha(4))
    plot(1:n,bigtna2.P,'k-')
    hold on;
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
  ylim([0 1])

        axes(ha(5))
    plot(1:n,bigtso2.Fe.*1000,'k-')
    hold on;
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('[Fe]')
    xlim([xl(1) xl(2)])
    ylim([0 0.3])

    axes(ha(6))
    plot(1:n,bigtna2.Fe.*1000,'k-')
    hold on;
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
    ylim([0 0.3])


      axes(ha(7))
    plot(1:n,bigtso2.muNS.*86400,'k-')
    hold on;
    plot(1:n,bigtso2.muPS.*86400,'b-')
    plot(1:n,bigtso2.muFeS.*86400,'r-')
    plot(1:n,bigtso2.muCS.*86400,'g-')
    plot(1:n,bigtso2.mumaxS.*86400,'m--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    l=legend({'N','P','Fe','C','max'},'numcolumns',5,...
        'EdgeColor','none','BackgroundAlpha',0);
    l.Position = [0.0747    0.71   0.4093    0.006];
    ylabel('\mu_{Cyano}')
    xlim([xl(1) xl(2)])
    ylim([0 2])
    % ylim([-0.2 2])


    axes(ha(8))
    plot(1:n,bigtna2.muNS.*86400,'k-')
    hold on;
    plot(1:n,bigtna2.muPS.*86400,'b-')
    plot(1:n,bigtna2.muFeS.*86400,'r-')
    plot(1:n,bigtna2.muCS.*86400,'g-')
    plot(1:n,bigtna2.mumaxS.*86400,'m--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
    ylim([0 2])


    axes(ha(9))
    plot(1:n,bigtso2.mub./365,'k-')
    hold on;
    plot(1:n,bigtso2.mu_f./365,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('\mu_B,\mu_F')
    xlim([xl(1) xl(2)])
    ylim([-0.2 0.8])

    axes(ha(10))
    plot(1:n,bigtna2.mub./365,'k-')
    hold on;
    plot(1:n,bigtna2.mu_f./365,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
    ylim([-0.2 0.8])

    axes(ha(11))
    plot(1:n,bigtso2.prodB,'k-')
    hold on;
   % plot(1:n,bigtso2.mortb,'r--')
    plot(1:n,bigtso2.mortb,'b--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('prod_B,mort_B')
    xlim([xl(1) xl(2)])
 ylim([0 4])
% ylim([0 1e-5])

    axes(ha(12))
    plot(1:n,bigtna2.prodB,'k-')
    hold on;
    plot(1:n,bigtna2.mortb,'r--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
  ylim([0 4])
 % ylim([0 1e-5])

    axes(ha(13))
    % -prod_p + remin_p + regen_p + tr.dop/taup;
    %plot(1:n,bigtso2.reminP+bigtso2.regenP+bigtso2.DOP-bigtso2.prodP,'k-')
    plot(1:n,bigtso2.B,'k-')
    hold on;
     plot(1:n,bigtso2.F,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('B, F')
    xlim([xl(1) xl(2)])
   ylim([0 0.03])

    axes(ha(14))
  %  plot(1:n,bigtna2.reminP+bigtna2.regenP+bigtna2.DOP-bigtna2.prodP,'k-')
      plot(1:n,bigtna2.B,'k-')
    hold on;
     plot(1:n,bigtna2.F,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
  ylim([0 0.03])


    axes(ha(15))
    plot(1:n,bigtso2.QnS,'k-')
    hold on;
    plot(1:n,bigtso2.QnL,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('Q_N')
    xlim([xl(1) xl(2)])
    ylim([0 0.2])

    axes(ha(16))
    plot(1:n,bigtna2.QnS,'k-')
    hold on;
    plot(1:n,bigtna2.QnL,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
    ylim([0 0.2])

    axes(ha(17))
    plot(1:n,bigtso2.QpS,'k-')
    hold on;
    plot(1:n,bigtso2.QpL,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('Q_P')
    xlim([xl(1) xl(2)])
   ylim([0 0.03])
   yline(QpmaxS,'r--')
   yline(QpmaxL,'r--')

    axes(ha(18))
    plot(1:n,bigtna2.QpS,'k-')
    hold on;
    plot(1:n,bigtna2.QpL,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
 ylim([0 0.03])
    yline(QpmaxS,'r--')
   yline(QpmaxL,'r--')

    axes(ha(19))
    plot(1:n,(bigtso2.VpS./bigtso2.QpS).*dpery,'k-')
    hold on;
    plot(1:n,(bigtso2.VpL./bigtso2.QpL).*dpery,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('\mu_{P}')
    xlim([xl(1) xl(2)])
 %  ylim([1 2])
 ylim([0 800])


    axes(ha(20))
    plot(1:n,(bigtna2.VpS./bigtna2.QpS).*dpery,'k-')
    hold on;
    plot(1:n,(bigtna2.VpL./bigtna2.QpL).*dpery,'k--')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
   ylim([0 800])

    axes(ha(21))
    plot(1:n,bigtso2.QnS./bigtso2.QpS,'k-')
    hold on;
    plot(1:n,bigtso2.QnL./bigtso2.QpL,'k--')
    plot(1:n,bigtso2.ro,'r-')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    ylabel('R_{NP}')
    xlim([xl(1) xl(2)])
   ylim([5 40])


    axes(ha(22))
    plot(1:n,bigtna2.QnS./bigtna2.QpS,'k-')
    hold on;
    plot(1:n,bigtna2.QnL./bigtna2.QpL,'k--')
    plot(1:n,bigtna2.ro,'r-')
    grid on; box on;
    set(gca,'xticklabel',{},'fontsize',8)
    xlim([xl(1) xl(2)])
    ylim([5 40])
end

    %%