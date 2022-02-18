% Bayesian Inversion, Level Set Inversion 1D
% Author: Kristina Backes
% 2022
close all; clear;
tic

%% Model

% Profil set up
AB = [1, 1.3, 1.8, 2.4, 3.2, 4.2, 5.6, 7.5, ...
        10, 13, 18, 24, 32, 42, 56, 75, ...
        100, 130, 180 240 320 420 560 750 ...
        1000 1300 1800 2400 3200 4200]; 
    
% Resisitivity
rho.back = 100; % Background
rho.ano = 1000; % Anomaly

% Set up model
Hom.rho = rho.back * ones(5,1);
Hom.rho(1) = rho.ano;
Hom.thk = [10 10 10 9];

Het.rho = rho.back * ones(5,1);
Het.rho(2) = rho.ano ;
Het.thk = [10 10 10 9];
len_split = length(Hom.rho);

% Add random noise
Het.rhoa = dcfwdf(Het.rho, Het.thk, AB);

% noise_mode: 
% 1: additive noise, 2: relative noise, 3: superposition 
noise_mode = 2;

switch noise_mode
    case 1
        % Additive noise
        noise = 3;
        Het.rhoa_n = Het.rhoa + noise .* randn(size(Het.rhoa));
        sigma = noise; % Standardabweichung
    case 2
        % Relative noise
        noise = 0.03;
        Het.rhoa_n = Het.rhoa .* (1 + noise * randn(size(Het.rhoa)) );
        sigma = noise; % Varianz
    case 3
        % Superposition relativ and additive
        noise = 0.03;
        noise_g = 0.003;
        noise_geo = noise_g * AB + zeros(size(AB));
        Het.rhoa_nn = Het.rhoa .* (1 + noise * randn(size(Het.rhoa)) );
        Het.rhoa_n = Het.rhoa_nn + noise_geo' .* randn(size(Het.rhoa));
        sigma = noise; % Varianz
end

% Number of observations
nn = length(Het.rhoa_n);

% Depth in meter
T = 40;

%% Define MCMC

% Number iterations
n = 35000; 

% Burn amount
burn = 8000;

% Standard deviation of normal proposal density or step size
s = 0.2;

%% KL-process
% Mesh points = center of 1m resolution
x = 0.5:1:(T+0.5);
K = length(x);

% Distance between mesh points
R = sqrt((x'-x).^2);

% Point by point covariance matrix
c = @(r, corrLen) exp(- r/corrLen);
C = c(R, 8);

% Eigensystem of covariance matrix
[V, Lambda] = eig(C);
L = V*sqrt(Lambda);

%% pCN - Algorithmus

% Random Walk
thk_curr = Hom.thk;
rhoa_curr = Hom.rho;
xi_curr = randn(K,1);

% Auxiliary variable
hlp_time = 1;
hlp = 0;
mittelA = zeros(1, n);
mittelwert_neu = zeros(K,1);
xi_prop = zeros(K,1);
post = zeros(1,n);

% Start of Random Walk
post(1) = 0;
target_u(:,1) = 0.5 + zeros(K,1); 
target_Fu(:,1) = 0.5 + zeros(K,1);
target_xi(:,1) = xi_curr;

% MCMC loop
for i = 2:n
        
    % Random numbers   
    xi_curr = target_xi(:,i-1);
    for k = 1:K
        xi_prop(k,1) = sqrt(1-s^2) * xi_curr(k,1) + s*randn(1,1);
    end
    
    u = L * xi_prop;

    % Tresholding binary vector
    Fu = (u >= 0);

    % Transform arrays for Het.rho and Het.thk 
    thk_prop = diff( [0; find( diff(Fu) )] )';
    rhoa_prob = rho.back + (rho.ano - rho.back) * Fu(cumsum(thk_prop));
    
    if length( unique(Fu) ) == 1
        A = 0;
        hlp = hlp + 1;
        
    else
        % Forward operation
        fwd_prop = dcfwdf(rhoa_prob, thk_prop, AB); % rhoa für proposed_x

        % Likelihood candidate
        switch noise_mode
            case 1
                % Additive noise
                log_like_prop = - (1/(2*sigma^2)) * ...
                    sum(((Het.rhoa_n)/1 - (fwd_prop)/1).^2);
            case 2
                % Relative noise
                log_like_prop = - sum(log(fwd_prop)) - (1/(2*sigma^2)) * ...
                    sum((Het.rhoa_n - fwd_prop).^2 ./ fwd_prop.^2);
            case 3
                % Superposition relativ and additive
                log_like_prop = - sum(log(fwd_prop * noise + noise_geo')) - ...
                    0.5 * sum((Het.rhoa_n - fwd_prop).^2 ./ ...
                    (fwd_prop.^2 * noise^2 + noise_geo'.^2));
        end
        
        % If previous step was accepted
        if hlp_time == 1  
            % Forward operation
            fwd_curr = dcfwdf(rhoa_curr, thk_curr, AB); % rhoa für current_x
            
            % Likelihood previous value i-1
            switch noise_mode
                case 1
                    % Additive noise
                    log_like_curr = - (1/(2*sigma^2)) * ...
                        sum(((Het.rhoa_n)/1 - (fwd_curr)/1).^2);
                case 2
                    % Relative noise
                    log_like_curr = - sum(log(fwd_curr)) - (1/(2*sigma^2)) * ...
                        sum((Het.rhoa_n - fwd_curr).^2 ./ fwd_curr.^2); 
                case 3
                    % Superposition relativ and additive
                    log_like_curr = - sum(log(fwd_curr * noise + noise_geo')) - ...
                        0.5 * sum((Het.rhoa_n - fwd_curr).^2 ./ ...
                        (fwd_curr.^2 * noise^2 + noise_geo'.^2));
            end 
        end
        
        % Probability of move
        A = exp(log_like_prop - log_like_curr);
    end
    
    % Generate uniformly distributed random variable
    uA = rand(1,1);
    
    % Acceptance query
    if uA <= min(A, 1)
        % Random walk
        thk_curr = thk_prop;
        rhoa_curr = rhoa_prob;
        target_u(:,i) = u; 
        target_Fu(:,i) = Fu;
        target_xi(:,i) = xi_prop;
        
        % Posterior density
        post(i) = exp(log_like_prop - 0.5 * sum(xi_prop)^2);
              
        % Results
        result.u(i-1,:) = u;
        result.Fu(i-1,:) = Fu;
        result.xi(i-1,:) = xi_prop;
        result.stat(i-1,1) = log_like_prop + log(post(i));
        
        % Auxiliary variable
        mittelA(i) = 1; 
        hlp_time = 1;
        
        % Online mean calculation
        % M_n = (n-1)/n * M_(n-1) + 1/n * x_n
        mittelwert_alt = mittelwert_neu;
        mittelwert_neu = (i-1)/i .* mittelwert_alt + 1/i .* Fu;
       
    else
        % Random walk       
        target_u(:,i) = target_u(:,i-1); 
        target_Fu(:,i) = target_Fu(:,i-1);
        target_xi(:,i) = target_xi(:,i-1);
        
        % Auxiliary variable
        mittelA(i) = 0; 
        hlp_time = 0;
        
        % Posterior density
        post(i) = post(i-1);
        
        % Online mean calculation
        % M_n = (n-1)/n * M_(n-1) + 1/n * x_n
        mittelwert_alt = mittelwert_neu;
        mittelwert_neu = (i-1)/i .* mittelwert_alt + 1/i .* target_Fu(:,i-1);        
    end
    
    As(i) = A;
end

Am = unique(As);
Am = Am';

% Acceptance probability of random walk
akwrs = mean(mittelA) * 100; 

%% Estimator 1 - MAP

% Sort for maximal posteriori density
[~, s_idx] = sort(result.stat(:, end), 'descend');
result.stat = result.stat(s_idx, :);
result.u = result.u(s_idx, :);
result.Fu = result.Fu(s_idx, :);

idx = find(result.stat(:,1), 1);
result.stat = result.stat(idx:end, :);
result.u = result.u(idx:end, :);
result.Fu = result.Fu(idx:end, :);

% Computation of the map_n most probable modell answers
map_n = 2;
modell_u = zeros(len_split, map_n);

for i = 1:map_n
    modell_u(:,i) = result.u(i,:);
end

% Statistics
modell_mean_u = mean(modell_u,2);
modell_mean_Fu = (modell_mean_u >= 0);

%% Estimator 2 - Mean of Markov Chain
% Statistics
modell_mean_u2 = mean(target_u(:,burn:end),2);
modell_mean_Fu2 = (modell_mean_u2 >= 0);

Fpm = mean(target_Fu(:,burn:end),2);

%% Plot

% Plot MAP
figure(1)
subplot(1,3,1)
imagesc(modell_mean_u(1:end-1))
title('A ---- MAP u')
colorbar

subplot(1,3,2)
imagesc(modell_mean_Fu(1:end-1))
title('A ---- MAP u --> Fu')
colorbar

% Plot Mean
figure(2)
subplot(1,3,1)
imagesc(modell_mean_u2(1:end-1))
title('B ---- Mean u')
colorbar

subplot(1,3,2)
imagesc(modell_mean_Fu2(1:end-1))
title('B ---- Mean u --> Fu')
colorbar

% Plot FPM
figure(3)
subplot(1,2,1)
imagesc(mittelwert_neu(1:end-1))
title('C1 ---- FPM - Mittelwert online')
colorbar

subplot(1,2,2)
imagesc(Fpm(1:end-1)) 
title('C2 ---- FPM') 
colorbar

%% Plot deep sounding comparison

Fu_hlp = (modell_mean_u >= 0);
thk_res = diff( [0; find( diff(Fu_hlp) )] )';
rhoa_res = rho.back + (rho.ano - rho.back) * Fu_hlp(cumsum(thk_res));
    
hlp_res = dcfwdf(rhoa_res, thk_res, AB);

figure(1)
subplot(1,3,3)
set(gca, 'Fontsize', 15)
loglog(AB, Het.rhoa_n, '--k', 'LineWidth', 1.2);
hold on
loglog(AB, hlp_res, '-b', 'LineWidth', 0.5);
title('A ---- Sondierungskurve MAP ---- A','FontWeight','bold','FontSize',12)
ylabel('\rho_a in \Omegam')
xlabel('AB/2 in m')
legend('Observed data','Result of MAP-estimator','location','northwest')

%% Mean u
Fu_hlp = (modell_mean_u2 >= 0);
thk_res = diff( [0; find( diff(Fu_hlp) )] )';
rhoa_res = rho.back + (rho.ano - rho.back) * Fu_hlp(cumsum(thk_res));
    
hlp_res = dcfwdf(rhoa_res, thk_res, AB);

figure(2)
subplot(1,3,3)
set(gca, 'Fontsize', 15)
loglog(AB, Het.rhoa_n, '--k', 'LineWidth', 1.2);
hold on
loglog(AB, hlp_res, '-b', 'LineWidth', 0.5);
title('B ---- Sondierungskurve von mean u','FontWeight','bold','FontSize',12)
ylabel('\rho_a in \Omegam')
xlabel('AB/2 in m')
legend('Observed data','Result of MAP-estimator','location','northwest')

%% Plot for MA

gg = rho.back* ones(size(Fu));
gg(11:21) = rho.ano;

% Tresholding binary vector
Fgg = (gg >= 500);

% Transform arrays for Het.rho and Het.thk 
thk_prop = diff( [0; find( diff(Fu) )] )';
rhoa_prob = rho.back + (rho.ano - rho.back) * Fgg(cumsum(thk_prop));
    
figure(4)

subplot(2,2,1)
imagesc(gg(1:end-1)) 
set(gca, 'Fontsize', 15)
title('Model') 
ylabel('Tiefe [m]')
setcolorbartitle('\rho [\Omegam]')

subplot(2,2,2)
imagesc(Fpm(1:end-1)) 
set(gca, 'Fontsize', 15)
ylabel('Tiefe [m]')
title('FPM') 
setcolorbartitle('%')

subplot(2,2,[3 4])
loglog(AB, Het.rhoa_n, 'o-r', 'LineWidth', 1.5);
hold on
loglog(AB, hlp_res, 'x-b', 'LineWidth', 1.5);
title('Sondierungskurve','FontWeight','bold')
set(gca, 'Fontsize', 15)
ylabel('\rho_a [\Omegam]')
xlabel('AB/2 [m]')
legend('Messwerte','Ergebnis MAP','location','northwest')
grid on
x0=100;
y0=1000;
width=720;
height=360;

%% Time display
et = toc;
etchar = sprintf('Elapsed time: %.2f sec. \n', et);
fprintf(etchar);

%% Helper

function setcolorbartitle(str)
    colorcet('D1')
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = str;
    set(colorTitleHandle ,'String',titleString);
end










