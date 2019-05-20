load('../data/cylinderDNS.mat','U','V')
m = size(U);
UALL = reshape(U, m(1)*m(2), m(3));
VALL = reshape(V, m(1)*m(2), m(3));
X_o = [reshape(U, m(1)*m(2), m(3)); reshape(V, m(1)*m(2), m(3))];
mkdir error_data
%% calculate errors
lambda = 1;
eta = [0:5:100]/100; % corruption fractions
eta1 = eta/2; % eta1 represents half the corruption percentage
eta2 = 1 - eta1; % eta2 accounts for the second half of the corruption
u = UALL;
v = VALL;

iters = 5; % number of points to converge to statistically significant results
 
for k = 1:length(eta1)
    
    parfor m = 1:iters
 
        
        % this is in a for loop to gather statistically significant results 
        % UALL and VALL are redefined as the original, uncorrupted data 
        UALL = u;
        VALL = v;
        
        
        rep = std(UALL(:))*10; % the value of the randomly chosen pixels
        
        x = rand(size(u)); % a randomly generated matrix of size u
        
        % replace the values in the UALL with +/- rep
        UALL(x<eta1(k)) = rep;
        UALL(x>eta2(k)) = -rep;
        
        % replace the same values in the VALL with +/- rep
        VALL(x<eta1(k)) = rep;
        VALL(x>eta2(k)) = -rep;
        
        %concatenate corrupted UALL and VALL
        X = [UALL; VALL];
        
       
        
        [n1,n2] = size(X);
        thresh = 1e-7;
        lam = lambda*1/sqrt(max(n1,n2));
        [L, S, count] = inexact_alm_rpca(X, lam, thresh, 1000);

        error(m) = norm(X_o-L, 'fro')/norm(X_o, 'fro');
        [X_u_temp, X_sigma, ~] = svd(X_o,'econ');
        X_u(:,:,m) = X_u_temp(:,1:20);
        
        [lowrank_u_temp,lowrank_sigma,~] = svd(L, 'econ');
        lowrank_u(:,:,m) = lowrank_u_temp(:,1:20);
        
        beta = min(size(L))/max(size(L));
        
        lowrank_sigma = diag(lowrank_sigma);
       
%         lowrank_sigma = lowrank_sigma(lowrank_sigma > 10^(-8));
%         [rank_L(m), tau] = getrank_beta(beta,lowrank_sigma);
        rank_L(m) = sum(lowrank_sigma)/sum(diag(X_sigma));
        
        [~, X_noise_sigma, ~] = svd(X,'econ');
        rank_L_noise(m) = sum(lowrank_sigma)/sum(diag(X_noise_sigma));
        %
        
    end
    filename = sprintf('error_data/error_density_%07.5f.mat', 2*eta1(k));
    
    save(filename, 'error', 'rank_L', 'rank_L_noise', 'lowrank_u', 'X_u')
    clearvars -except lambda d tol eta2 k m i j u v X_o eta1 iters
    
    
end


%% Gather all the error data
% initialize variables
error_all = [];
rank_all = [];
rank_noise_all = [];

for i = 1:length(eta)
    filename = sprintf('error_data/error_density_%07.5f.mat', eta(i));
    load(filename, 'error', 'rank_L', 'rank_L_noise')
    error_all = [error_all, error'];
    rank_all = [rank_all, rank_L'];
    rank_noise_all = [rank_noise_all, rank_L_noise'];
    %     for k = 1:5
    %     for j = 1:10
    %         X_all(j,k,i) = norm(X_u(:,j, k) - lowrank_u(:,j, k));
    %     end
    %     end
end



yyaxis right

plot(eta,median(error_all), '.', 'MarkerSize', 15); hold on
plot(eta,median(error_all), 'Linewidth', 1.2)
ylabel('Mean Relative Error')
ylim([-.05 1.05])

yyaxis left

plot(eta, mean(rank_all),'.', 'MarkerSize', 15); hold on
plot(eta, mean(rank_all), 'Linewidth', 1.2)
grid on

ylim([-.05 1.05])
ylabel('Relative Nuclear Norm') % left y-axis 
xlabel('% Corrupted Pixels')
