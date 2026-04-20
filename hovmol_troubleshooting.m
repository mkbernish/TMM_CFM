
%%
gr = grid;
figure('Position',[1 1 900 600])
[ha,pos]=tight_subplot(4,1,[0.08 0.016],[.12 .12],[.06 .1]);
axes(ha(1))
%clvl = linspace(0,5,12);
contourf(1:length(roH),gr.yt,roH,24,'linecolor','none')
shading interp;
ylim([-70 70])
xlim([1 60])
c1 = colorbar;
c1.Position = [0.91    0.68   0.01   0.2];
c1.Label.String = ['[NO_3^-] (\muM)'];
colormap(m_colmap('jet'))
clim([8 24])
set(gca,'fontname','times','fontsize',10)
ytickformat('degrees')
grid on; box on;

axes(ha(2))
contourf(1:length(roH),gr.yt,PH,24,'linecolor','none')
shading interp;
ylim([-70 70])
xlim([1 60])
c1 = colorbar;
c1.Position = [0.91    0.4   0.01   0.2];
c1.Label.String = ['[NO_3^-] (\muM)'];
colormap(m_colmap('jet'))
clim([0 1])
set(gca,'fontname','times','fontsize',10)
ytickformat('degrees')
grid on; box on;

axes(ha(3))
contourf(1:length(roH),gr.yt,BH,24,'linecolor','none')
shading interp;
ylim([-70 70])
xlim([1 60])
c1 = colorbar;
c1.Position = [0.91    0.12   0.01   0.2];
c1.Label.String = ['[NO_3^-] (\muM)'];
colormap(m_colmap('jet'))
clim([0 0.1])
set(gca,'fontname','times','fontsize',10)
ytickformat('degrees')
grid on; box on;

axes(ha(4))
contourf(1:length(roH),gr.yt,FeH,24,'linecolor','none')
shading interp;
ylim([-70 70])
xlim([1 60])
c1 = colorbar;
c1.Position = [0.91    0.12   0.01   0.2];
c1.Label.String = ['[NO_3^-] (\muM)'];
colormap(m_colmap('jet'))
clim([0 0.001])
set(gca,'fontname','times','fontsize',10)
ytickformat('degrees')
grid on; box on;


