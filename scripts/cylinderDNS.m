%% load and format data
load('../data/cylinderDNS.mat','U','V')
m = size(U);
X = [reshape(U, m(1)*m(2), m(3)); reshape(V, m(1)*m(2), m(3))];
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
Lu01 = reshape(L01(1:end/2, : ), m(1), m(2), m(3)); Lv01 = reshape(L01(1+end/2:end, : ), m(1), m(2), m(3));
Su01 = reshape(S01(1:end/2, : ), m(1), m(2), m(3)); Sv01 = reshape(S01(1+end/2:end, : ), m(1), m(2), m(3));

Lu02 = reshape(L02(1:end/2, : ), m(1), m(2), m(3)); Lv02 = reshape(L02(1+end/2:end, : ), m(1), m(2), m(3));
Su02 = reshape(S02(1:end/2, : ), m(1), m(2), m(3)); Sv02 = reshape(S02(1+end/2:end, : ), m(1), m(2), m(3));

Lu05 = reshape(L05(1:end/2, : ), m(1), m(2), m(3)); Lv05 = reshape(L05(1+end/2:end, : ), m(1), m(2), m(3));
Su05 = reshape(S05(1:end/2, : ), m(1), m(2), m(3)); Sv05 = reshape(S05(1+end/2:end, : ), m(1), m(2), m(3));

Lu1 = reshape(L1(1:end/2, : ), m(1), m(2), m(3)); Lv1 = reshape(L1(1+end/2:end, : ), m(1), m(2), m(3));
Su1 = reshape(S1(1:end/2, : ), m(1), m(2), m(3)); Sv1 = reshape(S1(1+end/2:end, : ), m(1), m(2), m(3));

Lu2 = reshape(L2(1:end/2, : ), m(1), m(2), m(3)); Lv2 = reshape(L2(1+end/2:end, : ), m(1), m(2), m(3));
Su2 = reshape(S2(1:end/2, : ), m(1), m(2), m(3)); Sv2 = reshape(S2(1+end/2:end, : ), m(1), m(2), m(3));

Lu5 = reshape(L5(1:end/2, : ), m(1), m(2), m(3)); Lv5 = reshape(L5(1+end/2:end, : ), m(1), m(2), m(3));
Su5 = reshape(S5(1:end/2, : ), m(1), m(2), m(3)); Sv5 = reshape(S5(1+end/2:end, : ), m(1), m(2), m(3));

Lu10 = reshape(L10(1:end/2, : ), m(1), m(2), m(3)); Lv10 = reshape(L10(1+end/2:end, : ), m(1), m(2), m(3));
Su10 = reshape(S10(1:end/2, : ), m(1), m(2), m(3)); Sv10 = reshape(S10(1+end/2:end, : ), m(1), m(2), m(3));

%% Plot results
load ../data/plotparams.mat cylinderDNScax cylinderDNScmap
rows = 2; cols = 7; i = 0;


figure(1);
for ind = 1:m(3)
    i = 0;
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu01(:,:,ind),Lv01(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 0.1')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu02(:,:,ind),Lv02(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 0.2')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu05(:,:,ind),Lv05(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 0.5')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu1(:,:,ind),Lv1(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 1')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu2(:,:,ind),Lv2(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    colormap(CC); title('L, \lambda = 2')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu5(:,:,ind),Lv5(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 5')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Lu10(:,:,ind),Lv10(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('L, \lambda = 10')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    % sparse start
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su01(:,:,ind),Sv01(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off; title('S, \lambda = 0.1')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su02(:,:,ind),Sv02(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 0.2')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su05(:,:,ind),Sv05(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 0.5')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su1(:,:,ind),Sv1(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 1')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su2(:,:,ind),Sv2(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 2')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su5(:,:,ind),Sv5(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 5')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    i = i+1;
    subplot(rows, cols, i)
    imagesc(flipud((curl(Su10(:,:,ind),Sv10(:,:,ind)))));
    shading flat; axis equal; axis tight; axis off;
    title('S, \lambda = 10')
    caxis(cylinderDNScax); colormap(cylinderDNScmap)
    
    drawnow
    
end
