%% load data and create X matrix
filename = []; %set your filename and path to your data formatted withu and v as velocities
load(filename,'U','V') % u and v are your velocities
m = size(U);
X = [reshape(U, m(1)*m(2), m(3)); reshape(V, m(1)*m(2), m(3))]; %reshape so columns are time and rows are space

%% add noise--for no noise eta = 0

eta = 0; % Set eta to your fraction of occluded vectors 
rep = std(X(:))*10;
x = rand(size(X(1:end/2,:)));
b = sort(x(:));
thresh = b(floor(.5*eta*numel(b)));

X(x<thresh) = rep;

x = rand(size(X(1:end/2,:)));
b = sort(x(:));
thresh = b(floor(.5*eta*numel(b)));

X(x<thresh) = -rep;


%% Run RPCA
lambda = 1; % choose your sparsity constant value (1 is often a good starting point)
tol = 1e-7; % set your tolerance
maxIter = 1000; % set your maximum number of iterations
[L, S, ~] = inexact_alm_rpca(X, lambda, tol, maxIter); %L is lowrank, S is sparse

%% Reshape your U and V low rank and sparse data

Lu = reshape(L(1:end/2, : ), m(1), m(2), m(3)); Lv = reshape(L(1+end/2:end, : ), m(1), m(2), m(3));
Su = reshape(S(1:end/2, : ), m(1), m(2), m(3)); Sv = reshape(S(1+end/2:end, : ), m(1), m(2), m(3));

%% Calculate turbulence Spectra

%for original data 
[k, E] = turbspec(U, V, W); % if you do not have a third velocity dimension, use zeroes(size(U))
%for low rank
[k_L, E_L] = turbspec(Lu, Lv, Lw); % if you do not have a third velocity dimension, use zeroes(size(Lu))
