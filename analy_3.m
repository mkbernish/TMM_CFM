addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\"))
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
%reg  = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_cfdi__cmub_varFe_10_reg10_arr_zprod');
reg = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_cfdi__cmub_varFe_15_reg11_test');
par = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_cfdi__cmub_varFe_10_reg10_arr_par_zprod');
temp = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_cfdi__cmub_varFe_10_reg10_arr_temp_zprod');
rfe = load('tr_all_DOM_Dfree_rdi8_diazFe_1_explic_cfdi__cmub_varFe_10_reg10_arr_zprod_rfe2p.mat');

%%
plt.biomass = 1;
if plt.biomass
    figure('Position',[800 200 700 600])
    [ha,pos]=tight_subplot(3,1,[0.04 0.04],[.1 .1],[.1 .05]);
    axes(ha(1))
    plot(reg.biomass_np,'k-o','markerfacecolor','k')
    grid on; box on;
    hold on;
    plot(par.biomass_np,'r-o','markerfacecolor','r')
    plot(temp.biomass_np,'b-o','markerfacecolor','b')
    plot(rfe.biomass_np,'g-o','markerfacecolor','g')
    % ylim([15 20])
    ylabel('Biomass N:P')
    set(gca,'xticklabel',{})
    % xlabel('Year')

    axes(ha(2))
    plot(reg.N_inv_surf./reg.P_inv_surf,'k-o','markerfacecolor','k')
    grid on; box on;
    hold on;
    plot(par.N_inv_surf./par.P_inv_surf,'r-o','markerfacecolor','r')
    plot(temp.N_inv_surf./temp.P_inv_surf,'b-o','markerfacecolor','b')
    plot(rfe.N_inv_surf./rfe.P_inv_surf,'g-o','markerfacecolor','g')
    %  ylim([11.2 12.6])
    set(gca,'xticklabel',{})
    ylabel('Surface \SigmaN: \SigmaP')

    axes(ha(3))
    plot(reg.N_inv./reg.P_inv,'k-o','markerfacecolor','k')
    grid on; box on;
    hold on;
    plot(par.N_inv./par.P_inv,'r-o','markerfacecolor','r')
    plot(temp.N_inv./temp.P_inv,'b-o','markerfacecolor','b')
    plot(rfe.N_inv./rfe.P_inv,'g-o','markerfacecolor','g')
    ylim([14.2 14.3])
    ylabel('\SigmaN: \SigmaP')
    xlabel('Year')
    legend({'Baseline','Increased PAR','Increased temp','R_{Fe:P}'},'numcolumns',3)
end
%%
plt.np =1;
if plt.np
    fs = 11; zz = 1;
    ww = 0*M3d+NaN; ww(iocn) =temp.ro-reg.ro;
    cmap = flipud(getPyPlot_cMap('RdBu', 256,'pyCmd','py'));
    figure('Position',[100 400 950 500])
    [ha,pos]=tight_subplot(2,2,[0.04 0.04],[.06 .06],[.04 .06]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(-3,3,24);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)),xl,'linecolor','none')
    hold on;shading interp;
    clim([-3 3])
    % clim([-0.01 0.01])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    %  cb1.Position = [ 0.12    0.77   0.32   0.024];
    cb1.Label.String = '\Delta R_O ';
    %  cb1.Ticks = [8, 12, 16, 20, 24, 28, 32];
    %
    axes(ha(2))
    ww = 0*M3d+NaN; ww(iocn) =temp.outL.fracBiosynth-reg.outL.fracBiosynth;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(-0.2,0.2,20);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)),xl,'linecolor','none')
    hold on;shading interp;
    clim([-0.2 0.2])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);
    % cb2.Position = [ 0.59    0.77   0.32   0.024];
    cb2.Label.String = '\Delta Biosynth fraction';
    %  cb2.Ticks = [0, 0.25,0.5,0.75,1];

    axes(ha(3))
    ww = 0*M3d+NaN; ww(iocn) =temp.tr.QnL-reg.tr.QnL;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(-0.04,0.04,20);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)),xl,'linecolor','none')
    hold on;shading interp;
    clim([-0.04 0.04])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);
    % cb2.Position = [ 0.59    0.77   0.32   0.024];
    cb2.Label.String = '\Delta N:C';
    %  cb2.Ticks = [0, 0.25,0.5,0.75,1];

    axes(ha(4))
    ww = 0*M3d+NaN; ww(iocn) =temp.tr.QpL-reg.tr.QpL;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(-0.004,0.004,20);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)),xl,'linecolor','none')
    hold on;shading interp;
    clim([-0.004 0.004])
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'yticklabel',{},'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb2=colorbar('fontsize',fs);
    % cb2.Position = [ 0.59    0.77   0.32   0.024];
    cb2.Label.String = '\Delta P:C';
    %  cb2.Ticks = [0, 0.25,0.5,0.75,1];

end

%% limitation
plt.limt = 1;
if plt.limt
    fs = 11; zz = 1;
    ww = 0*M3d+NaN; ww(iocn) =reg.outL.limType;
    cmap = flipud(getPyPlot_cMap('Pastel1', 4,'pyCmd','py'));
    figure('Position',[100 400 900 500])
    [ha,pos]=tight_subplot(1,2,[0.04 0.04],[.06 .06],[.04 .06]);
    axes(ha(1))
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(1,4,4);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,zz),3,'omitnan'),45,2)),xl,'linecolor','none')
    hold on;shading interp;
    clim([1 4])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    cb1.Ticks = linspace(1.5,3.5,4);
    cb1.TickLabels = {'N','P','Fe','C'};


    axes(ha(2))
    ww = 0*M3d+NaN; ww(iocn) =temp.outL.limType;
    m_proj('miller','lat',[-80 80],...
        'lon',[-180 180]);
    xl = linspace(1,4,4);
    m_contourf(grid.xt-180,grid.yt,...
        (circshift(mean(ww(:,:,zz),3,'omitnan'),45,2)),xl,'linecolor','none')
    hold on;shading interp;
    clim([1 4])
    colormap(cmap(:,1:3))
    m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
    m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
    cb1=colorbar('fontsize',fs);
    cb1.Ticks = linspace(1.5,3.5,4);
    cb1.TickLabels = {'N','P','Fe','C'};
end

%%
fs = 11; zz = 1;
ww = 0*M3d+NaN; ww(iocn) =reg.tr.Fe.*1000;
%cmap = flipud(getPyPlot_cMap('RdBu', 256,'pyCmd','py'));
figure('Position',[100 400 700 400])
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.06 .06],[.04 .06]);
axes(ha(1))
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
% xl = linspace(-3,3,24);
xl = linspace(0,2,24);
m_contourf(grid.xt-180,grid.yt,...
    (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)),xl,'linecolor','none')
hold on;shading interp;
clim([0 2])
% clim([-0.01 0.01])
%  colormap(cmap(:,1:3))
cmocean('deep')
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
%  cb1.Position = [ 0.12    0.77   0.32   0.024];
cb1.Label.String = '\Delta R_O ';
%  cb1.Ticks = [8, 12, 16, 20, 24, 28, 32];
























%% determining compensation depth?
ww = 0*M3d+NaN; ww(iocn) =obs.I;
perc1 = ww(:,:,1).*0.01;
msk1 = 0*M3d;
cc = 0;
for i = 1 : size(msk1,3)
    cc = cc + 1;
    ok = ww(:,:,i)<=perc1;
    msk1(:,:,cc) = ok;
end
[~,idif2] = setdiff(iocn,find(msk1==0)); % no productivity subregion

%% Fe : P
fs = 12;
cmap = flipud(getPyPlot_cMap('RdYlBu', 256,'pyCmd','py'));
ww = 0*M3d+NaN; ww(iocn) =tr.QfeS./tr.QpS;
fep_d = mean(ww(:,:,1:2),[2 3],'omitnan');
figure('Position',[100 400 700 400])
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.12 .12],[.1 .03]);
axes(ha(1))
plot(grid.yt,fep_d,'k-')
yline(rfe2p,'k--')
xlim([-75 75])
grid on; box on;

%% global Fe: C
fs = 11; zz = 2;
M3d = reg.M3d; iocn = reg.iocn; gr = reg.grid;
ww = 0*M3d+NaN; ww(iocn) =reg.outL.FeC*1000000;
femax =1.44e-4.*1000000;
cmap = flipud(getPyPlot_cMap('RdYlBu', 256,'pyCmd','py'));
figure('Position',[100 200 900 650])
[ha,pos]=tight_subplot(3,2,[0.04 0.04],[.06 .06],[.04 .06]);
axes(ha(1))
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_pcolor(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)))
hold on;shading interp;
m_contour(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,2),3,'omitnan'),45,2)),[femax femax],'linecolor','g')
clim([0 120])
colormap(cmap(:,1:3))
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
cb1.Label.String = 'Fe:C (mol mol^{-1}) ';

axes(ha(2))
ww = 0*M3d+NaN; ww(iocn) =reg.outL.limType;
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_pcolor(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)))
hold on;shading interp;
colormap(cmap(:,1:3))
clim([1 4])
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
cb1.Label.String = 'Fe:C (mol mol^{-1}) ';


axes(ha(3))
ww = 0*M3d+NaN; ww(iocn) =reg.tr.QfeL./reg.tr.QpL;
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_pcolor(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)))
hold on;shading interp;
m_contour(gr.xt-180,gr.yt,...
    (circshift(mean(Pmat(:,:,2),3,'omitnan'),45,2)),[reg.rfe2p reg.rfe2p],'linecolor','g')
clim([0 0.04])
colormap(cmap(:,1:3))
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
cb1.Label.String = 'Fe:P (mol mol^{-1}) ';


axes(ha(4))
ww = 0*M3d+NaN; ww(iocn) =reg.tr.Fe.*1000;
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_pcolor(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)))
hold on;shading interp;
m_contour(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,1),3,'omitnan'),45,2)),[0 0],'linecolor','g')
clim([0 2])
colormap(cmap(:,1:3))
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
cb1.Label.String = '[Fe] nM ';

axes(ha(5))
ww = 0*M3d+NaN; ww(iocn) =reg.ro;
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_pcolor(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,1:2),3,'omitnan'),45,2)))
hold on;shading interp;
m_contour(gr.xt-180,gr.yt,...
    (circshift(mean(Pmat(:,:,2),3,'omitnan'),45,2)),[0 0],'linecolor','g')
clim([8 28])
m_scatter(gr.xt(8)-180,gr.yt(26),60,'r','filled')
colormap(cmap(:,1:3))
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
cb1.Label.String = 'Q_P^{diat} (mol P mol C^{-1}) ';

axes(ha(6))
ww = 0*M3d+NaN; ww(iocn) =reg.prod_b;
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_pcolor(gr.xt-180,gr.yt,...
    (circshift(mean(ww(:,:,1:zz),3,'omitnan'),45,2)))
hold on;shading interp;
m_contour(gr.xt-180,gr.yt,...
    (circshift(mean(Pmat(:,:,2),3,'omitnan'),45,2)),[0 0],'linecolor','g')
clim([0 4])
colormap(cmap(:,1:3))
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xtick',[-120, -60, 0, 60, 120],'xticklabel',{},'yticklabel',{});
cb1=colorbar('fontsize',fs);
cb1.Label.String = 'prod_b ';

%% effect of A_Pho_Fe on Qfe and muFe
wwn = 0*M3d+NaN; wwn(iocn) =reg.tr.N;
wwp = 0*M3d+NaN; wwp(iocn) =reg.tr.P;
afe = linspace(60/1000000,3.8e-3,200);
      [outS] =  alloc2(80,32,0.2,1,...
           0.1/1000,PmaxS,knS,kpS,kfeS,vmaxnS,vmaxpS,vmaxfeS,...
            QpmaxS,afe,EaS);

      figure('Position',[100 200 600 550])
[ha,pos]=tight_subplot(2,1,[0.04 0.04],[.06 .06],[.1 .1]);
axes(ha(1))
plot(afe,outS.mu_Fe,'k-')
hold on;
yline(outS.mu_N,'k--')
yline(outS.mu_P,'b--')
ylim([0 4e-5])

axes(ha(2))
plot(afe,outS.FeC,'k-')
ylim([0 4e-5])
