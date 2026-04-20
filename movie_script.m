%% making movies for Ro, N, P, and fraction

%% N
gr = grid;
Nmat = 0*M3d+NaN; Nmat(iocn) = tr.N;
Pmat = 0*M3d+NaN; Pmat(iocn) = tr.P;
romat = 0*M3d+NaN; romat(iocn) = ro;
fsmat = 0*M3d+NaN; fsmat(iocn) = tr.B;
mudimat = 0*M3d+NaN; mudimat(iocn) = tr.fL;
fs = 10; sz = 40;
figure('Position',[1 1 800 500],'visible','off')
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.04 .04],[.06 .1]);
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_contourf(gr.xt-180,gr.yt,circshift...
    (Nmat(:,:,1),45,2),44,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
clim([0 32])
cmocean('balance',44)
cb1=colorbar('fontname','times new roman','fontsize',fs);
cb1.Position = [ 0.91    0.24   0.015   0.5];
cb1.Label.String = ['[N]'];
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'Nmap_' num2str(i) '.jpeg']);
% "C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_coupling\Ninv_Fe_v2\cfdi__omub"
pth = ['C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_ss\Ninv_Fe_v2\' mstr '\N\'];
saveas(gcf,[pth fig_name])
close;

%% P
figure('Position',[1 1 800 500],'visible','off')
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.04 .04],[.06 .1]);
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_contourf(gr.xt-180,gr.yt,circshift...
    (Pmat(:,:,1),45,2),44,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
clim([0 2])
cmocean('balance',44)
cb1=colorbar('fontname','times new roman','fontsize',fs);
cb1.Position = [ 0.91    0.24   0.015   0.5];
cb1.Label.String = ['[P]'];
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'Pmap_' num2str(i) '.jpeg']);
pth = ['C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_ss\Ninv_Fe_v2\' mstr '\P\'];
saveas(gcf,[pth fig_name])
close;
%% ro
figure('Position',[1 1 800 500],'visible','off')
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.04 .04],[.06 .1]);
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_contourf(gr.xt-180,gr.yt,circshift...
    (romat(:,:,1),45,2),44,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
clim([0 40])
cmocean('balance',44)
cb1=colorbar('fontname','times new roman','fontsize',fs);
cb1.Position = [ 0.91    0.24   0.015   0.5];
cb1.Label.String = ['N:P'];
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'romap_' num2str(i) '.jpeg']);
pth = ['C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_ss\Ninv_Fe_v2\' mstr '\ro\'];
saveas(gcf,[pth fig_name])
close;



%% prodN
figure('Position',[1 1 800 500],'visible','off')
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.04 .04],[.06 .1]);
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_contourf(gr.xt-180,gr.yt,circshift...
    (fsmat(:,:,1),45,2),44,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
clim([0 0.1])
cmocean('balance',44)
cb1=colorbar('fontname','times new roman','fontsize',fs);
cb1.Position = [ 0.91    0.24   0.015   0.5];
cb1.Label.String = ['f_{di}'];
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'fsmap_' num2str(i) '.jpeg']);
pth = ['C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_ss\Ninv_Fe_v2\' mstr '\prodn\'];
saveas(gcf,[pth fig_name])
close;

%% mu_di
figure('Position',[1 1 800 500],'visible','off')
[ha,pos]=tight_subplot(1,1,[0.04 0.04],[.04 .04],[.06 .1]);
m_proj('miller','lat',[-80 80],...
    'lon',[-180 180]);
m_contourf(gr.xt-180,gr.yt,circshift...
    (mudimat(:,:,1),45,2),44,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[1 1 1],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
clim([0 1])
cmocean('balance',44)
cb1=colorbar('fontname','times new roman','fontsize',fs);
cb1.Position = [ 0.91    0.24   0.015   0.5];
cb1.Label.String = ['\mu_{di}'];
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'mudimap_' num2str(i) '.jpeg']);
pth = ['C:\Users\mathp\OneDrive\Documents\Papers\paper4\tmm_ss\Ninv_Fe_v2\' mstr '\mu_di\'];
saveas(gcf,[pth fig_name])
close;
