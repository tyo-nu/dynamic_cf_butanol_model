function CI = CI_matrix(M, p, dim)
% Calculate CI for each dim of M - 
% dim==1 is CI for each column (default), dim==2 is CI for each row - 

if ~exist('dim','var')
    dim = 1;
end

if dim == 2
    M = M';
end

n_vectors = size(M,2);
CI = zeros(2,n_vectors);

for i=1:n_vectors
    single_CI = CI_vector(M(:,i), p);
    CI(:,i) = single_CI(:);
end

if dim == 2
    CI = CI';
end
end


function CI = CI_vector(x, p)

% Taken from Adam Danz - 
% https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval

CI_alt = @(x,p) std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * ...
    tinv(abs([0,1]-(1-p)/2),sum(~isnan(x(:)))-1) + mean(x(:),'omitnan'); 

% Function to get percentil (z-test) CI
CI_z = @(x,p)prctile(x,abs([0,100]-(100-p*100)/2));

% Now taken from https://www.mathworks.com/help/stats/tinv.html

% Compute the sample mean, standard error, and degrees of freedom.
xbar = mean(x);
n = sum(~isnan(x(:)));
se = std(x)/sqrt(n);
nu = n - 1;
% Find the upper and lower confidence bounds for the 95% confidence interval.
alpha = 1 - p;
alphaLo = alpha/2;
alphaHi = 1 - alpha/2;
% Compute the critical values for the confidence bounds.
crit = tinv([alphaLo alphaHi], nu);
% Determine the confidence interval for the population mean.
CI = xbar + crit*se;

% Instead use bootstrapping - doesn't assume normal distribution
% [CI_mean_std, bootstat] = bootci(1e2, {@(x)[mean(x) std(x)], x}, ...
%     'Alpha', alpha,...
%     'Type', 'bca',...
%     'NBootStd',10);
% 
% CI = CI_mean_std(:, 1);

% Can also just do raw percentiles

CI = prctile(x, (100.*(1 - [alphaHi alphaLo])) );

end