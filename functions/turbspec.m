function [k, E] = turbspec(uvel, vvel, wvel)
% uvel = U.u; vvel = U.v; wvel = U.w;

%% Find dimensions
[ny, nx, m] = size(uvel);
n = nx*ny;          % Snapshot size

uvel = reshape(uvel, n, m);  % Reshape so columns are snapshots, rows are time series
vvel = reshape(vvel, n, m);
wvel = reshape(wvel, n, m);
mask = any(isnan(uvel), 2) & any(isnan(vvel), 2) & any(isnan(wvel), 2);  % Get rid of any invalid points
uvel(isnan(uvel)) = 0; vvel(isnan(vvel)) = 0; wvel(isnan(wvel)) = 0;

% Mean flow
u_mean = mean(uvel, 2);  v_mean = mean(vvel, 2); w_mean = mean(wvel, 2);

%% Time-averaged turbulence KE
E = zeros(ny, nx);
for i=1:m
    u = uvel(:, i) - u_mean;
    u_hat = fft2(reshape(u, ny, nx));    
    v = vvel(:, i) - v_mean;
    v_hat = fft2(reshape(v, ny, nx));
    w = wvel(:, i) - w_mean;
    w_hat = fft2(reshape(w, ny, nx));
    E = E + 0.5*(u_hat.*conj(u_hat) + v_hat.*conj(v_hat) + w_hat.*conj(w_hat));
end
E = E/m;

%% Wavenumber calculations
Lx = nx; Ly = ny;
kx=(2*pi/Lx)*[1:(nx/2) (-nx/2):-1]'; kx(1)=10^(-6);
kx = fftshift(kx);
ky=(2*pi/Ly)*[1:(ny/2) (-ny/2):-1]'; ky(1)=10^(-6);
ky = fftshift(ky);
E = fftshift(E);  % This is now the energy spectrum as a function of the 2D wavenumber (kx, ky)

[KX, KY] = meshgrid(kx, ky);

%% Polar integral to get energy spectrum  
% !!!! may need to retool this if in 3D
[~, KR] = cart2pol(KX, KY);   % Radial values of k over the mesh

k_range = linspace(0, (max(ky) - min(ky))/2, 100);  % Discretize k as a radial variable
dk = k_range(2)- k_range(1);

%% Mask KR and bin PSD to do a polar integral
E_k = zeros(length(k_range), 1);
for i=1:length(k_range)
    k = k_range(i);
    k_mask = (KR >= k) & (KR < k+dk);  % Identify locations in this bin
    E_disk = abs(E(k_mask));
    dtheta = 2*pi/sum(sum(k_mask));         % Approximate difference between points in the bin
    E_k(i) = 2*k^2*sum(E_disk(:))*dtheta;   % Riemann integral approximation
end

k = k_range;  E = E_k;  % Rename for output