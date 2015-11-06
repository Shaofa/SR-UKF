function [X, S] = srufk_update(Xpre, Spre, Y, h, Rs, h_param, sigmax, Wm, Wc)
    L = size(Xpre, 1);
    O = size(Y, 1);
    
    % pass sigmax through measure module
    sigmay = zeros(O, 2*L+1);        % 3x13
    Ypre = zeros(O, 1);
    for i = 1:2*L+1
        sigmay(:,i) = feval(h, sigmax(:,i));
        Ypre = Ypre + Wm(i) * sigmay(:, i);
    end
    % measure predict
    %Ypre = sigmay*Wm;
    
   % measure covirance predict
%     tmp = sqrt(Wc(2)) * (sigmay(:,2:2*L+1) - repmat(Ypre,1,2*L));
    for i= 2 : size(sigmay, 2)
        tmp(:, i-1) = sqrt(Wc(2)) * (sigmay(:, i) - Ypre);
    end
    [tmp, Sy] = qr([tmp Rs]', 0);
    Sy = cholupdate( Sy, sqrt(Wc(1)) * (sigmay(:,1) - Ypre) );
    
    Pxy = zeros(L, O);
    for i = 1 : 2*L+1
        Pxy = Pxy + Wc(i) * (sigmax(:,i) - Xpre) * (sigmay(:,i) - Ypre)';
    end
    
    K = Pxy/Sy/Sy';
    X = Xpre + K*(Y-Ypre);
    U = K*Sy';
    for i = 1:O
        Spre = cholupdate(Spre, U(:,i), '-');
    end
    S = Spre;
end