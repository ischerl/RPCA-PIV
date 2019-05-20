%% load and format data
load('../data/channelDNS.mat','U','V','W'); 
m = size(U);
X = [reshape(U, m(1)*m(2), m(3)); reshape(V, m(1)*m(2), m(3)); reshape(W, m(1)*m(2), m(3))];
%% add salt & pepper
eta = 0; % to add noise change eta to fraction of occluded points. eta = 0.2 for 20%
if eta ~=0
    rep = std(X(:))*10;
    x = rand(size(X(1:end/2,:)));
    b = sort(x(:));
    thresh = b(floor(.5*eta*numel(b)));
    
    X(x<thresh) = rep;
    
    x = rand(size(X(1:end/2,:)));
    b = sort(x(:));
    thresh = b(floor(.5*eta*numel(b)));
    
    X(x<thresh) = -rep;
end
%% take RPCA
[L01, S01, ~] = inexact_alm_rpca(X, 0.1, 1e-7, 1000);
[L02, S02, ~] = inexact_alm_rpca(X, 0.2, 1e-7, 1000);
[L05, S05, ~] = inexact_alm_rpca(X, 0.5, 1e-7, 1000);
[L1, S1, ~] = inexact_alm_rpca(X, 1, 1e-7, 1000);
[L2, S2, ~] = inexact_alm_rpca(X, 2, 1e-7, 1000);
[L5, S5, ~] = inexact_alm_rpca(X, 5, 1e-7, 1000);
[L10, S10, ~] = inexact_alm_rpca(X, 10, 1e-7, 1000);

%% reshape to get U & V for each
Lu01 = reshape(L01(1:end/3, : ), m(1), m(2), m(3)); Lv01 = reshape(L01(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw01 = reshape(L01(1+2*end/3:end, : ), m(1), m(2), m(3));
Su01 = reshape(S01(1:end/3, : ), m(1), m(2), m(3)); Sv01 = reshape(S01(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw01 = reshape(S01(1+2*end/3:end, : ), m(1), m(2), m(3));

Lu02 = reshape(L02(1:end/3, : ), m(1), m(2), m(3)); Lv02 = reshape(L02(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw02 = reshape(L02(1+2*end/3:end, : ), m(1), m(2), m(3));
Su02 = reshape(S02(1:end/3, : ), m(1), m(2), m(3)); Sv02 = reshape(S02(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw02 = reshape(S02(1+2*end/3:end, : ), m(1), m(2), m(3));

Lu05 = reshape(L05(1:end/3, : ), m(1), m(2), m(3)); Lv05 = reshape(L05(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw05 = reshape(L05(1+2*end/3:end, : ), m(1), m(2), m(3));
Su05 = reshape(S05(1:end/3, : ), m(1), m(2), m(3)); Sv05 = reshape(S05(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw05 = reshape(S05(1+2*end/3:end, : ), m(1), m(2), m(3));

Lu1 = reshape(L1(1:end/3, : ), m(1), m(2), m(3)); Lv1 = reshape(L1(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw1 = reshape(L1(1+2*end/3:end, : ), m(1), m(2), m(3));
Su1 = reshape(S1(1:end/3, : ), m(1), m(2), m(3)); Sv1 = reshape(S1(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw1 = reshape(S1(1+2*end/3:end, : ), m(1), m(2), m(3));

Lu2 = reshape(L2(1:end/3, : ), m(1), m(2), m(3)); Lv2 = reshape(L2(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw2 = reshape(L2(1+2*end/3:end, : ), m(1), m(2), m(3));
Su2 = reshape(S2(1:end/3, : ), m(1), m(2), m(3)); Sv2 = reshape(S2(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw2 = reshape(S2(1+2*end/3:end, : ), m(1), m(2), m(3));

Lu5 = reshape(L5(1:end/3, : ), m(1), m(2), m(3)); Lv5 = reshape(L5(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw5 = reshape(L5(1+2*end/3:end, : ), m(1), m(2), m(3));
Su5 = reshape(S5(1:end/3, : ), m(1), m(2), m(3)); Sv5 = reshape(S5(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw5 = reshape(S5(1+2*end/3:end, : ), m(1), m(2), m(3));

Lu10 = reshape(L10(1:end/3, : ), m(1), m(2), m(3)); Lv10 = reshape(L10(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Lw10 = reshape(L10(1+2*end/3:end, : ), m(1), m(2), m(3));
Su10 = reshape(S10(1:end/3, : ), m(1), m(2), m(3)); Sv10 = reshape(S10(1+end/3:2*end/3, : ), m(1), m(2), m(3)); Sw10 = reshape(S10(1+2*end/3:end, : ), m(1), m(2), m(3));

%% Plot results
load ../data/plotparams.mat channelDNScmap channelDNScax

rows = 2; cols = 7; i = 0;

figure(1);
for ind = 1%:m(3) % to play all as a video, delete first %
    i = 0;
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu01(:,:,ind),Lv01(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 0.1')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu02(:,:,ind),Lv02(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 0.2')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu05(:,:,ind),Lv05(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 0.5')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu1(:,:,ind),Lv1(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 1')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu2(:,:,ind),Lv2(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 2')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu5(:,:,ind),Lv5(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 5')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Lu10(:,:,ind),Lv10(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 10')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    % sparse start
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su01(:,:,ind),Sv01(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off; title('S, \lambda = 0.1')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su02(:,:,ind),Sv02(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 0.2')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su05(:,:,ind),Sv05(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 0.5')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su1(:,:,ind),Sv1(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 1')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su2(:,:,ind),Sv2(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 2')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su5(:,:,ind),Sv5(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 5')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(((curl(Su10(:,:,ind),Sv10(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 10')
    caxis(channelDNScax); colormap(channelDNScmap)
    
    drawnow
    
end

%% calculate tke

[k0, E0] = turbspec(U, V, W);
[k01, E01] = turbspec(Lu01, Lv01, Lw01);
[k02, E02] = turbspec(Lu02, Lv02, Lw02);
[k05, E05] = turbspec(Lu05, Lv05, Lw05);
[k1, E1] = turbspec(Lu1, Lv1, Lw1);
[k2, E2] = turbspec(Lu2, Lv2, Lw2);
[k5, E5] = turbspec(Lu5, Lv5, Lw5);
[k10, E10] = turbspec(Lu10, Lv10, Lw10);

%% plot TKE
figure;
red = [1 0 0];
orange = [255 128 0]/255;
green = [0 204 0]/255;
blue = [0 0 255]/255;
purple = [153 0 153]/255;
grey = [128,128,128]/255;

loglog(k0, E0,'k-', 'linewidth', 5); hold on
loglog(k05, E05,'color',purple, 'linewidth', 2);
loglog(k1, E1,'color',blue, 'linewidth', 2);
loglog(k2, E2,'color',green, 'linewidth', 2);
loglog(k5, E5,'color',orange, 'linewidth', 2);
loglog(k10, E10,'color',red, 'linewidth', 2);
titlename = sprintf('%d percent', eta*100);
title(titlename)
ylim([100-15 10^6])
xlim([min(k0) max(k0)+.2])

legend('original', '\lambda = 0.5','\lambda = 1','\lambda = 2','\lambda = 5','\lambda = 10', 'Location','southwest')
grid on
xlabel('k')
ylabel('E')
