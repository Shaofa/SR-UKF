function main(flag)
    clc;
    clear;
    close all;
    
    if(nargin < 1)
        file = '../data/filter_in.txt';
    else
        file = '../data/simulation.txt';
    end
    
%% get Measurements
    [raIn, xyIn, time, n] = getMeasure(file);
    raOut = zeros(size(raIn));
    xyOut = zeros(size(xyIn));
    
%% initialize srukf structure
    srukf.X = xyIn(:,1);
    srukf.L = size(srukf.X, 1);
    srukf.O = 3;
    %srukf.S = eye(srukf.L);
    srukf.S = chol(diag([0.6, 0.9, 0.1, 0.4, 0.3, 0.07]));
    srukf.Qs = chol(diag([0.2, 0.02, 0.002, 0.2, 0.02, 0.002]));
    srukf.Rs = chol(diag([0.5, 0.00087, 0.5]));
    srukf.alpha = 0.1;
    srukf.beta = 2.0;
    srukf.kappa = 0;
    srukf.time = 0.0;
    srukf.errCnt = 0;
    srukf.f = @srukf_f_func;
    srukf.h = @srukf_h_func;
    [srukf.Wm,srukf.Wc] = srukf_ut_weights(srukf.L,srukf.alpha, srukf.beta, srukf.kappa);
    
    %% filtering
    hasError = 0;
    for i = 1 : n
        [Xpre, Spre, sigmax] = srukf_predict( ....
            srukf.X,...
            srukf.S,...
            srukf.f,...
            srukf.Qs,...
            time(i) - srukf.time,...
            srukf.Wm,...
            srukf.Wc,...
            srukf.alpha,...
            srukf.beta,...
            srukf.kappa);
        raPre = srukf.h(Xpre);
        error = abs(raPre - raIn(:,i));
        if( error(1) > 2 || error(2) > 3*pi/180)
             srukf.errCnt = srukf.errCnt + 1;
             hasError = 1;
        else
            hasError = 0;
            srukf.errCnt = 0;
        end
        if(srukf.errCnt > 4)
            srukf.errCnt = 0;
            srukf.X = xyIn(:,i);
        end
        if(~hasError)
            [srukf.X, srukf.S] = srufk_update(...
                Xpre, ...
                Spre,...
                raIn(:,i),...
                srukf.h,...
                srukf.Rs,...
                [],...
                sigmax,...
                srukf.Wm,...
                srukf.Wc);
        end
        srukf.time = time(i);
        xyOut(:, i) = srukf.X;
        raOut(:, i) = srukf.h(srukf.X);
    end

    %% plot
    figure;
    plot(xyIn(1,:), xyIn(4,:), '*b');
    hold on;
    plot(xyOut(1,:), xyOut(4,:),'-r', 'linewidth',2);
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    xlabel('x/m');
    ylabel('y/m');
    legend('Measurement', 'Smoothed');
    title('SR-UKF');
    
    figure;
    subplot(2,1,1);
    plot(raIn(1,:),'-b', 'LineWidth', 1);
    hold on;
    plot(raOut(1,:),'-r',  'LineWidth', 2);
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    ylabel('Range/m');
    legend('Measurement', 'Smoothed');
    title('SR-UKF');
    
    subplot(2,1,2);
    plot(90-raIn(2,:)*180/pi, '-b', 'LineWidth', 1);
    hold on;
    plot(90-raOut(2,:)*180/pi, '-r', 'LineWidth', 2);
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    legend('Measurement', 'Smoothed');
    ylabel('Angle/degree');
    xlabel('time/0.048 sec');
    
end

%% get measurements
function [raIn, xyIn, time, n] = getMeasure(path)
    raw = dlmread(path)';
    n = size(raw, 2);
    raIn = zeros(3, n);
    time = raw(1, :);
    raIn(1, :) = raw(2,:);
    raIn(2, :) = (90 - raw(3,:)) .* pi ./ 180;
    raIn(3, :) = raw(4,:);
    xyIn(1, :) = raIn(1, :) .* cos(raIn(2, :));
    xyIn(2, :) = raIn(3, :) .* cos(raIn(2, :));
    xyIn(3, :) = 0;
    xyIn(4, :) = raIn(1, :) .* sin(raIn(2, :));
    xyIn(5, :) = raIn(3, :) .* sin(raIn(2, :));
    xyIn(6, :) = 0;
end

%% Dynamic model function
function [ x_n ] = srukf_f_func(x, dt)
    F=[ 1,  dt, dt*dt/2.,   0,  0,  0;          ...
        0,  1,  dt,         0,  0,  0;         ...
        0,  0,  1,          0,  0,  0;         ...
        0,  0,  0,          1,  dt, dt*dt/2.;  ...
        0,  0,  0,          0,  1,  dt;        ...
        0,  0,  0,          0,  0,  1];
    x_n = F*x;
end

%% Measure model function
function [ y ] = srukf_h_func(x)
    y(1,:) = sqrt(x(1,:).^2+x(4,:).^2);
    y(2,:) = atan2(x(4,:),x(1,:));
    y(3,:) = (x(1,:).*x(2,:)+x(4,:).*x(5,:))./y(1,:);
end

%% UT weights
function [WM, WC] = srukf_ut_weights(n, alpha, beta, kappa)
    lambda = alpha^2 * (n + kappa) - n; 
    WM = zeros(2*n+1,1);
    WC = zeros(2*n+1,1);
    for j=1:2*n+1
        if j==1
            WM(j) = lambda / (n + lambda);
            WC(j) = lambda / (n + lambda) + (1 - alpha^2 + beta);
        else
            WM(j) = 1 / (2 * (n + lambda));
            WC(j) = WM(j);
        end
    end
end