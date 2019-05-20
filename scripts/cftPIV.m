%% Crossflow turbine PIV

%% cross-correlated

% load data and reshape into X
load('../data/cftPIV_cc.mat','U','V', 'phase')
U_cc = U; V_cc = V;
m = size(U_cc);
X = double([reshape(U_cc, m(1)*m(2), m(3)); reshape(V_cc, m(1)*m(2), m(3))]);
rep = 10*std(U_cc(:));

k = find(X == 0);

r = randi([0 1], length(k), 1);

X(k(r == 0)) = rep;
X(k(r == 1)) = -rep;
% calculate RPCA
[L, S, ~] = inexact_alm_rpca(X, 1.6, 1e-7, 1000);

% Reshape low rank (L) and sparse (S) into U and V
Lu_cc = reshape(L(1:end/2, : ), m(1), m(2), m(3)); Lv_cc = reshape(L(1+end/2:end, : ), m(1), m(2), m(3));
Su_cc = reshape(S(1:end/2, : ), m(1), m(2), m(3)); Sv_cc = reshape(S(1+end/2:end, : ), m(1), m(2), m(3));


%% vector validated

% load data and reshape into X
load('../data/cftPIV_mf.mat','U','V')
U_mf = U; V_mf = V;
m = size(U_mf);
X = double([reshape(U_mf, m(1)*m(2), m(3)); reshape(V_mf, m(1)*m(2), m(3))]);
rep = 10*std(U_mf(:));

k = find(X == 0);

r = randi([0 1], length(k), 1);

X(k(r == 0)) = rep;
X(k(r == 1)) = -rep;
% calculate RPCA
[L, S, ~] = inexact_alm_rpca(X, 1.6, 1e-7, 1000);

% Reshape low rank (L) and sparse (S) into U and V
Lu_mf = reshape(L(1:end/2, : ), m(1), m(2), m(3)); Lv_mf = reshape(L(1+end/2:end, : ), m(1), m(2), m(3));
Su_mf = reshape(S(1:end/2, : ), m(1), m(2), m(3)); Sv_mf = reshape(S(1+end/2:end, : ), m(1), m(2), m(3));

%% interpolated

% load data and reshape into X
load('../data/cftPIV_int.mat','U','V');
U_int = U; V_int = V;
m = size(U_int);
X = double([reshape(U_int, m(1)*m(2), m(3)); reshape(V_int, m(1)*m(2), m(3))]);

% calculate RPCA
[L, S, ~] = inexact_alm_rpca(X, 1.6, 1e-7, 1000);

% Reshape low rank (L) and sparse (S) into U and V
Lu_int = reshape(L(1:end/2, : ), m(1), m(2), m(3)); Lv_int = reshape(L(1+end/2:end, : ), m(1), m(2), m(3));
Su_int = reshape(S(1:end/2, : ), m(1), m(2), m(3)); Sv_int = reshape(S(1+end/2:end, : ), m(1), m(2), m(3));

%% Calculate statistics
ph = linspace(0, pi, 33);

% add phase to data .mat files

for i = 1:length(ph)-1
   
   [ind, ~] = find((phase > ph(i)) & (phase <= ph(i+1)));
   
   tp_L_mf = Lu_mf(:,:,ind);
   tp_X_i = U_int(:,:,ind);
   
   pmed_L_mf(:,:,i) = nanmedian(tp_L_mf,3);
   pmed_X_i(:,:,i) = nanmedian(tp_X_i, 3);

   pstd_X_i(:,:,i) = nanstd(tp_X_i,0,3);
   pstd_L_mf(:,:,i) = nanstd(tp_L_mf,0,3);
end

index = [];
for i = 1:length(phase)
    index_temp = find(phase(i) >= ph);
    index = [index, max(index_temp)];
    
end



%% Plot RPCA results and statistics
load ../data/plotparams cftPIVintcmap cftPIVstdcax cftPIVstdcmap cftPIVcax cftPIVcmap
U_cc(U_cc==0) = NaN;
U_mf(U_mf==0) = NaN;

c = [-25, 25];
rows = 2;
cols = 5;

figure(1)
for ind = 1:200
    
    i = 1;
    ax1 = subplot(rows, cols,i);
    imagesc(fliplr(U_cc(:,:,ind)))
    caxis(ax1, cftPIVcax)
    colormap(ax1, cftPIVcmap)
    shading flat; axis equal; axis tight; axis off; 
    title('Cross-Correlated')
    
    i = i+1;
    ax2 = subplot(rows, cols,i);
    imagesc(fliplr(U_mf(:,:,ind)))
    shading flat; axis equal; axis tight; axis off;
    caxis(ax2, cftPIVcax); 
    colormap(ax2, cftPIVcmap); 
    title('Vector-Validated')
    
    
    i = i+1;
    ax3 = subplot(rows, cols,i);
    imagesc(fliplr(U_int(:,:,ind)));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax3, cftPIVcax); 
    colormap(ax3,cftPIVintcmap); 
    title('Interpolated')
    
    i = i+1;
    ax4 = subplot(rows, cols,i);
    imagesc(fliplr(pmed_X_i(:,:,index(ind))));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax4, cftPIVcax); 
    colormap(ax4, cftPIVcmap); 
    title('Interpolated, phase median')
    
    i = i+1;
    ax5 = subplot(rows, cols,i);
    imagesc(fliplr(pstd_X_i(:,:,index(ind))));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax5, cftPIVstdcax); 
    colormap(ax5, cftPIVstdcmap); 
    title('Interpolated, standard deviation')
    
    i = 1+i;
    ax6 = subplot(rows, cols,i);
    imagesc(fliplr(Lu_cc(:,:,ind)));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax6, cftPIVcax); 
    colormap(ax6, cftPIVcmap); 
    title('Cross-Correlated, L')
    
    i = i+1;
    ax7 = subplot(rows, cols,i);
    imagesc(fliplr(Lu_mf(:,:,ind)));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax7, cftPIVcax); 
    colormap(ax7, cftPIVcmap); 
    title('Vector-Validated, L')
    
    i = i+1;
    ax8 = subplot(rows, cols,i);
    imagesc(fliplr(Lu_int(:,:,ind)));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax8, cftPIVcax); 
    colormap(ax8, cftPIVcmap); 
    title('Interpolated, L')
    
    i = i+1;
    ax9 = subplot(rows, cols,i);
    imagesc(fliplr(pmed_L_mf(:,:,index(ind))));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax9, cftPIVcax); 
    colormap(ax9, cftPIVintcmap); 
    title('Vector Validated, L, phase median')
    
    i = i+1;
    ax10 = subplot(rows, cols,i);
    imagesc(fliplr(pstd_L_mf(:,:,index(ind))));
    shading flat; axis equal; axis tight; axis off;
    caxis(ax10, cftPIVstdcax); 
    colormap(ax10, cftPIVstdcmap); 
    title('Vector Validated, L, standard deviation')
    
    drawnow

end




