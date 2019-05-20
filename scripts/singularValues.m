 %% load data & take SVD

load('../data/channelDNS.mat','U','V','W')
m = size(U);
X = [reshape(U, m(1)*m(2), m(3)); reshape(V, m(1)*m(2), m(3)); reshape(W, m(1)*m(2), m(3))];
[~,channel_s,~] = svd(X, 'econ');

load('../data/cylinderPIV.mat','U','V')
m = size(U);
X = [reshape(U,m(1)*m(2), m(3)); reshape(V,m(1)*m(2), m(3))];
[~,cyl_piv_s,~] = svd(X, 'econ');

load('../data/cylinderDNS.mat','U','V')
m = size(U);
X = [reshape(U,m(1)*m(2), m(3)); reshape(V,m(1)*m(2), m(3))];
[~,cyl_sim_s,~] = svd(X, 'econ');

load('../data/cftPIV_cc.mat','U','V')
m = size(U);
X = [reshape(U,m(1)*m(2), m(3)); reshape(V,m(1)*m(2), m(3))];
[~,turb_piv_s,~] = svd(X, 'econ');

%% Plot Singular Values
ms = 10;
figure(1)
loglog(diag(channel_s)/sum(diag(channel_s)),'.', 'markersize', ms); hold on
loglog(diag(cyl_piv_s)/sum(diag(cyl_piv_s)),'.', 'markersize', ms)
loglog(diag(cyl_sim_s)/sum(diag(cyl_sim_s)),'.', 'markersize', ms)
loglog(diag(turb_piv_s)/sum(diag(turb_piv_s)),'.', 'markersize', ms)
legend('Channel Flow','Cylinder PIV','Cylinder DNS','Turbine PIV');
grid on


%% Plot Cumulative Energy
figure(2)
semilogx(cumsum(diag(channel_s)/sum(diag(channel_s))),'.', 'markersize', ms); hold on
semilogx(cumsum(diag(cyl_piv_s)/sum(diag(cyl_piv_s))),'.', 'markersize', ms)
semilogx(cumsum(diag(cyl_sim_s)/sum(diag(cyl_sim_s))),'.', 'markersize', ms)
semilogx(cumsum(diag(turb_piv_s)/sum(diag(turb_piv_s))),'.', 'markersize', ms)
ylim([0 1])
legend('Channel Flow','Cylinder PIV','Cylinder Simulation','Turbine PIV'); grid on
grid on