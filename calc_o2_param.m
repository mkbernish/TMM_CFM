% calculate o2 gas exchange parameters

clear all;
close all;
sperd=86400; dpery = 365.25; spery = sperd*dpery; % time conversions
    
%% paths
setpaths;
mod_mat = '90x180x24'; 
base_path = ['/work/tweber/MITgcm/Ninv_Fe/'];
model_path = ['/work/tweber/MITgcm/' mod_mat];
glacial = 0;

%% get transport matrix & grid
load([model_path '/Matrices/grid_matrices'],'M3d','grid','MSKS'); 
iocn = find(M3d(:)==1); 
m = length(iocn);
isurf = find(M3d(:,:,1)==1);
i3d = M3d*0+NaN; i3d(iocn) = 1:m;
[j1,isub,j2] = intersect(iocn,isurf);
[j1,idif] = setdiff(iocn,isurf);
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d;
AREA = grid.DXT3d.*grid.DYT3d;
dv = VOL(iocn);
da = AREA(iocn);

%% get data
load([model_path '/Data/WOA05_monthly']);
load([model_path '/Data/ocmip_gasex']);
temp_mon_glacial = temp_mon - 2;

%% calculate o2sat in mmol/m^3

if glacial;
    for l = 1:12;
        o2sat_star_mon(:,:,l) = o2sato(temp_mon_glacial(:,:,1,l),salt_mon(:,:,1,l));
    end

else
    for l = 1:12;
        o2sat_star_mon(:,:,l) = o2sato(temp_mon(:,:,1,l),salt_mon(:,:,1,l));
    end

end

o2sat_mon = p_mon.*o2sat_star_mon;

%% calculate sc_o2 (dimensionless)
for l = 1:12;
    Sc_o2_mon(:,:,l) = sc_o2(temp_mon(:,:,1,l));
end

%% calculate kw_o2

% in m/s
kw_o2_mon = (1-fice_mon).*xkw_mon.*sqrt(660./Sc_o2_mon);

% convert to m/yr
kw_o2_mon = kw_o2_mon*spery;

%% calculate annual averages

o2sat = squeeze(nanmean(o2sat_mon,3));
kw_o2 = squeeze(nanmean(kw_o2_mon,3));

%% save
if glacial; save([model_path '/Data/o2_gasex_glacial'],'kw_o2*','','o2sat*'); 
else save([model_path '/Data/o2_gasex'],'kw_o2*','','o2sat*'); end
