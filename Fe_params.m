%% Testing Fe parameters
% look at parameters for solving the Fe system

clear all;
close all;
sperd=86400; dpery = 365; spery = sperd*dpery; % time conversions

%% paths
setpaths;
mod_mat = '45x90x24'; % if 45x90x45, transport A4 gives better results than A5;
base_path = ['/work/tweber/MITgcm/Ninv_Fe/'];
model_path = ['/work/tweber/MITgcm/' mod_mat];

load([model_path '/Matrices/grid_matrices'],'M3d','grid','A4','MSKS'); 
A = A4; clear A4; 
m = size(A,1);
iocn = find(M3d(:)==1); 
isurface = find(M3d(:,:,1)==1);
i3d = M3d*0+NaN; i3d(iocn) = 1:m;
[j1,isurf,j2] = intersect(iocn,isurface);
[j1,idif] = setdiff(iocn,isurf);
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
AREA = grid.DXT3d.*grid.DYT3d;
dv = VOL(iocn);
da = AREA(iocn);

%% get data
load(fullfile(model_path, '/Data/WOA05_annual'));
load(fullfile(model_path, '/Data/Dust_1yr'));
dust = dust_mat(iocn);
dust4lim = repmat(dust_mat(:,:,1),[1 1 length(grid.zt)]);
d4lim = dust4lim(iocn);

%% Johnson (low ligand, high K)

% concentrations
FeT = [0:0.01:2];
LT = 0.6;

% equilibrium constant 
K = 1.2e13;

% solve
[Fe_free,L_free,FeL] = solve_Fe(FeT,LT,K);

% subplot(211);
% plot(FeT,Fe_free,'-b'); hold on;
% plot(FeT,L_free,'-r'); 
% legend('Location','North','free Fe','free Lligand');
% axis([0 2 0 1.5]);


%% Parekh / Follows (higher ligand, lower K)

% concentrations
FeT = [0:0.01:2];
LT = 1;

% equilibrium constant 
K = 2e11;

% solve
[Fe_free,L_free,FeL] = solve_Fe(FeT,LT,K);

% subplot(212);
% plot(FeT,Fe_free,'-b'); hold on;
% plot(FeT,L_free,'-r'); 
% legend('Location','North','free Fe','free Lligand');
% axis([0 2 0 1.5]);

%% Fe input rate from dust
fFe = 0.035; % fraction of iron in dust
Fe_in = fFe*(dust(isurf)'*da(isurf))/56; % mol/yr

%% calculate Fe loss
% follows params
LT = 1/1000; % mmol/m3
ksc = 1.1e-3*365;
K = 2e5;

[Fe_free,L_free,FeL] = solve_Fe(obs.Fe,LT,K);
Fe_scav = ksc*Fe_free; % mmol/m3/yr
Fe_out = (Fe_scav'*dv)/1000; % mol/yr

%% calculate bioavailable fraction
aFe = Fe_out/Fe_in;




