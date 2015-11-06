 function [Xpre, Spre, sigmax] = srukf_predict(X, S, f, Qs, f_param, Wm, Wc, alpha, beta, kappa)   
    L = size(X, 1);
    theta = sqrt(alpha^2 * (L + kappa));
    sigma = [zeros(size(X)), S, -S];
    sigma = repmat(X, 1, size(sigma,2)) + sigma*theta;
    %sigma = ut_sigmas(X, S, theta);
    
    % pass sigma through measure module
    % and state predict
    sigmax = zeros(size(sigma));            % 6x13
    Xpre = zeros(L, 1);
    for i = 1:size(sigma, 2)
        sigmax(:,i) = feval(f, sigma(:,i), f_param);
        Xpre = Xpre + Wm(i)* sigmax(:, i);
    end
  
    % covirance predict
    for i = 2 : size(sigmax, 2)
        tmp(:, i-1) = sqrt(Wc(2)) * (sigmax(:, i) - Xpre);
    end
    %tmp = sqrt(Wc(2)) * (sigmax(:,2:2*L+1) - repmat(Xpre,1,2*L));
    [tmp, Spre] = qr([tmp Qs]', 0);
    Spre = cholupdate(Spre,Wc(1)*(sigmax(:,1)-Xpre));
 end
 
%% sigma = [M M .. M] + sqrt(theta)[0 S S] 
function sigma = ut_sigmas(M, S, theta)
    sigma = [zeros(size(M)) S -S];
    sigma = repmat(M, 1, size(sigma,2)) + sqrt(theta)*sigma;
end