% Bayesian Inversion
% Author: Kristina Backes
% 5-layer system
% 2022

% With this script you can invert a 5-layer system with fixed thickness.
% dcfwf.m is required.

close all; clear; clc;

% rng(1, "twister");

%% Set up model 
tic

% profil
AB = [1, 1.3, 1.8, 2.4, 3.2, 4.2, 5.6, 7.5, ...
        10, 13, 18, 24, 32, 42, 56, 75, ...
        100, 130, 180 240 320 420 560 750];

% set up model
Hom.rho = [500 500 500 500 500];
Hom.thk = [1 5 5 10];
Het.rho = [100 50 50 1500 2000];
Het.thk = [1 5 5 10];
len_split = length(Hom.rho);

% add random noise
Het.rhoa = dcfwdf(Het.rho, Het.thk, AB);

% noise type: 
% 1: Additiver Noise, 2: Multipilakter noise, 3: Ueberlagerung
noise_type = 2;

switch noise_type
    case 1
        % Additiver Noise
        noise = 3;
        Het.rhoa_n = Het.rhoa + noise .* randn(size(Het.rhoa));
        sigma = noise; % Standardabweichung
    case 2
        % Multipilakter noise
        noise = 0.03;
        Het.rhoa_n = Het.rhoa .* (1 + noise * randn(size(Het.rhoa)) );
        sigma = noise; % Varianz
    case 3
        % Ueberlagerung multiplikativ und additiv
        noise = 0.03;
        noise_g = 0.003;
        noise_geo = noise_g * AB + zeros(size(AB));
        Het.rhoa_nn = Het.rhoa .* (1 + noise * randn(size(Het.rhoa)) );
        Het.rhoa_n = Het.rhoa_nn + noise_geo' .* randn(size(Het.rhoa));
        sigma = noise; % Varianz
end

nn = length(Het.rhoa_n);

%% Define MCMC

% number iterations
n = 39000; 

% burn = burnin amount
burn = 3000;

% standard deviation of normal proposal density or step size
s = 0.01; 

% Random walk
target = zeros(len_split, n);
target(:,1) = Hom.rho;
result.layer = zeros(n-1, len_split);
result.stat = zeros(n-1, 2);

% Auxiliary variable
rhoa_bestmodell = zeros(nn, 20);
modell = zeros(len_split, 20);
mittelA = 0;
prio_hist = zeros(1, n);

%% Metropolis-Hastings-Algorithmus

% Acceptance query
hlp_time = 1;

% Loop random walk
for i = 2:n
    
    % Assignement of previous position in MCMC chain
    current_log_rho = log(target(:,i-1));
    
    % Disturb every modell parameter independently 
    % Generating new candidate for random walk
    for k = 1:len_split 
        proposed_log_rho(k) = current_log_rho(k) + s*randn(1,1);
    end
         
    % Prior
    prior_prop = normpdf(proposed_log_rho, 5.5, 0.9);
    prior_curr = normpdf(current_log_rho, 5.5, 0.9);
     
    % Forward operation candidate
    fwd_prop = dcfwdf(exp(proposed_log_rho), Hom.thk, AB);

    % Likelihood candidate
    switch noise_type
        case 1
            % Additive Log-Like
            Log_Like_prop = - (1/(2*sigma^2)) * ...
                sum(((Het.rhoa_n)/1 - (fwd_prop)/1).^2);
        case 2
            % Multiplikative Log-Like
            Log_Like_prop = - sum(log(fwd_prop)) - (1/(2*sigma^2)) * ...
                sum((Het.rhoa_n - fwd_prop).^2 ./ fwd_prop.^2);
        case 3
            % Ueberlagerung
            Log_Like_prop = - sum(log(fwd_prop * noise + noise_geo')) - ...
                0.5 * sum((Het.rhoa_n - fwd_prop).^2 ./ ...
                (fwd_prop.^2 * noise^2 + noise_geo'.^2));
    end
    
    % Forward operation and likelihood, if previous step was accepted
    if hlp_time == 1  
        % Forward operation
        fwd_curr = dcfwdf(exp(current_log_rho), Hom.thk, AB);
        
        % Likelihood
        switch noise_type
            case 1
                % Additive Log-Like
                Log_Like_curr = - (1/(2*sigma^2)) * ...
                    sum(((Het.rhoa_n)/1 - (fwd_curr)/1).^2);
            case 2
                % Multiplikative Log-Like
                Log_Like_curr = - sum(log(fwd_curr)) - (1/(2*sigma^2)) * ...
                    sum((Het.rhoa_n - fwd_curr).^2 ./ fwd_curr.^2); 
            case 3
                % Ueberlagerung
                Log_Like_curr = - sum(log(fwd_curr * noise + noise_geo')) - ...
                    0.5 * sum((Het.rhoa_n - fwd_curr).^2 ./ ...
                    (fwd_curr.^2 * noise^2 + noise_geo'.^2));
        end 
    end
    
    % Posterior density
    post_dens = sum(log(prior_prop)) + Log_Like_prop;

    % Probability of move
    A = exp(post_dens - sum(log(prior_curr)) - Log_Like_curr); 
    
    % Generate uniformly distributed random variable
    u = rand(1,1);
    
    % Acceptance query
    if u <= min(A, 1)
        % Random walk
        target(:,i) = exp(proposed_log_rho); % accept proposal as new
        
        % Results
        result.layer(i-1,:) = exp(proposed_log_rho);
        result.stat(i-1,1) = Log_Like_prop;
        result.stat(i-1,2) = post_dens;
        
        % Auxiliary variable
        mittelA = mittelA + 1; % akzeptiert neuen Move
        hlp_time = 1;
    else
        % Random walk
        target(:,i) = exp(current_log_rho); % set old to be the new
        
        % Auxiliary variable
        hlp_time = 0;
    end
end

% Acceptance probability of random walk
akwrs = (mittelA/n) * 100; 

%% Priori-distribution 

for i = 1:n
    prio_hist(i) = lognrnd(5.5, 0.9); % mean = 244,7; sigma = 2,46
end

%% Estimator 1 - MAP

% Sort for maximal posteriori density
[~, s_idx] = sort(result.stat(:, end), 'descend');
result.stat = result.stat(s_idx, :); 
result.layer = result.layer(s_idx, :); 

idx = find(result.stat(:,1), 1);
result.stat = result.stat(idx:end, :);
result.layer = result.layer(idx:end, :);

% Computation of the 20 most probable modell answers
for i = 1:20
    rhoa_bestmodell(:,i) = dcfwdf(result.layer(i,:), Het.thk, AB);
    modell(:,i) = result.layer(i,:);
end

% Statistics
modell_mean = mean(modell,2);
modell_var = var(modell,0,2);
modell_std = std(modell,0,2);
modell_stat = [modell_mean modell_std];

% Quantil
rho_quant = quantile(modell,[0.025 0.25 0.50 0.75 0.975],2);

%% Estimator 2 - Mean of Markov Chain
% Statistics
modell_mean2 = mean(target(:,burn:end),2);
modell_var2 = var(target(:,burn:end),0,2);
modell_std2 = std(target(:,burn:end),0,2);
modell_stat2 = [modell_mean2 modell_std2];

% Quantil
rho_quant2 = quantile(target(:,burn:end),[0.025 0.25 0.50 0.75 0.975],2);

% Model answer
rhoa_bestmodell2 = dcfwdf(modell_mean2', Het.thk, AB);

%% Plot Random Walk
maxxi = max(target, [], 'all');
minni = min(target, [], 'all');

figure(1)
for m = 1:len_split
    subplot(2,3,m)
    plot(1:n, target(m,:), 'ko'); hold on
    plot(1:n, ones(1,n)*Het.rho(m), 'LineWidth', 2.5);
    plot([burn, burn], [minni, maxxi], '-b', 'LineWidth', 1.5);
    set(gca, 'Fontsize', 13)
    title("Schicht " + m + " (" + Het.rho(m) + " \Omegam)")
    set(gca, 'YScale', 'log')
    xlabel('Index')
    if m == 1 || m == 4
        ylabel('Zustand \rho [\Omegam]')
    end
    xlim([0,n])
    ylim([40, maxxi])
%     ylim([min(target(m,:)), max(target(m,:))])
    if m == 1
        legend('Zustand', 'Wahrer Wert', 'burn in')
    end
    grid on
end

%% Plot deep sounding comparison

figure(2)
loglog(AB, Het.rhoa_n, 'o-r', 'LineWidth', 1.5);
hold on
loglog(AB, rhoa_bestmodell2, 'x-b', 'LineWidth', 1);
for p = 1:20
    loglog(AB, rhoa_bestmodell(:,p), 'x-k', 'LineWidth', 1);
end
set(gca, 'Fontsize', 15)
title('Sondierungskurve','FontWeight','bold')
ylabel('\rho_a in \Omegam')
xlabel('AB/2 in m')
legend('Messwerte','Mittelwert MC','MAP Ergebnisse','location','northwest')
grid on
x0=100;
y0=1000;
width=720;
height=360;
set(gcf,'position',[x0,y0,width,height])

%% Comparison Priori & Posteriori
% Rho: 100 50 50 1500 2000

for k = 1:len_split
    [f(k,:), xi(k,:)] = ksdensity(target(k,:), (1:2:4000));
end

prior_histks = ksdensity(prio_hist, (1:2:4000));

mu = 5.5;
sig = 0.9;

m = exp(mu + sig^2/2);
v = exp(2*mu + sig^2) * (exp(sig^2)-1);
sig = sqrt(v);

figure(3)
for i = 1:5
    subplot(2,3,i)
    set(gca, 'Fontsize', 16)
    hold on
    plot(xi(1,:), f(i,:), 'LineWidth', 2.5)
    plot(xi(1,:), prior_histks, 'LineWidth', 2.5)
    set(gca, 'XScale', 'log')
    title("Schicht " +i+ " (" + Het.rho(i) + "\Omegam)")
    xlim([10 4000])
    xticks(logspace(1,4,4))
    xlabel('\rho [\Omegam]')
    legend("\mu_0","\pi", 'location', 'east' ,'Fontsize',16)
%     legend("\mu_0: \mu: " + round(m,1) + ", \sigma: " + round(sig,1),...
%         "\pi: \mu: " +...
%         round(modell_stat2(i,1),1) + ", \sigma: " + ...
%         round(modell_stat2(i,2),1) , 'location', 'northeast' ,'Fontsize',14)
     grid on
     ax = gca;
     ax.YAxis.Exponent = 0;
end

%% Plot grey zones 
% Rho: 100 50 50 1500 2000
% q: [0.025 0.25 0.50 0.75 0.975]
% [x y w h]

% Layer thickness
d = [1 5 5 10 9];
d_len = sum(d);

% Maximal expansion Plot
max_hlp = max(max(target(:)), 2100); 

figure(4)
set(gca, 'Fontsize', 16)
hold on

% Layer 1
rectangle('Position',[rho_quant2(1,1) 0 abs(minus(rho_quant2(1,1),rho_quant2(1,5))) 1], ...
    'FaceColor',[0.85 0.85 0.85])
rectangle('Position',[rho_quant2(1,2) 0 abs(minus(rho_quant2(1,2),rho_quant2(1,4))) 1], ...
    'FaceColor',[0.4 0.4 0.4])
plot(ones(1, length(0:1))*Het.rho(1), 0:1, 'r-', 'LineWidth', 1.5);
plot(50:100, ones(1, length(50:100)), 'r-', 'LineWidth', 1.5);

% Layer 2
rectangle('Position',[rho_quant2(2,1) 1 abs(minus(rho_quant2(2,1),rho_quant2(2,5))) 5], ...
    'FaceColor',[0.85 0.85 0.85])
rectangle('Position',[rho_quant2(2,2) 1 abs(minus(rho_quant2(2,2),rho_quant2(2,4))) 5], ...
    'FaceColor',[0.4 0.4 0.4])
plot(ones(1, length(1:6))*Het.rho(2), 1:6, 'r-', 'LineWidth', 1.5);

% Layer 3
rectangle('Position',[rho_quant2(3,1) 6 abs(minus(rho_quant2(3,1),rho_quant2(3,5))) 5], ...
    'FaceColor',[0.85 0.85 0.85])
rectangle('Position',[rho_quant2(3,2) 6 abs(minus(rho_quant2(3,2),rho_quant2(3,4))) 5], ...
    'FaceColor',[0.4 0.4 0.4])
plot(ones(1, length(6:11))*Het.rho(3), 6:11, 'r-', 'LineWidth', 1.5);
plot(50:1500, 11*ones(1, length(50:1500)), 'r-', 'LineWidth', 1.5);

% Layer 4
rectangle('Position',[rho_quant2(4,1) 11 abs(minus(rho_quant2(4,1),rho_quant2(4,5))) 10], ...
    'FaceColor',[0.85 0.85 0.85])
rectangle('Position',[rho_quant2(4,2) 11 abs(minus(rho_quant2(4,2),rho_quant2(4,4))) 10], ...
    'FaceColor',[0.4 0.4 0.4])
plot(ones(1, length(11:21))*Het.rho(4), 11:21, 'r-', 'LineWidth', 1.5);
plot(1500:2000, 21*ones(1, length(1500:2000)), 'r-', 'LineWidth', 1.5);
plot(50:1500, 11*ones(1, length(50:1500)), 'r-', 'LineWidth', 1.5);

% Layer 5
rectangle('Position',[rho_quant2(5,1) 21 abs(minus(rho_quant2(5,1),rho_quant2(5,5))) 9], ...
    'FaceColor',[0.85 0.85 0.85])
rectangle('Position',[rho_quant2(5,2) 21 abs(minus(rho_quant2(5,2),rho_quant2(5,4))) 9], ...
    'FaceColor',[0.4 0.4 0.4])
plot(ones(1, length(21:30))*Het.rho(5), 21:30, 'r-', 'LineWidth', 1.5);

set(gca,'XAxisLocation','top')
set(gca,'YDir','reverse');
set(gca, 'XScale', 'log')
xlim([0 max_hlp])
ylim([0 d_len])
xlabel('\rho [\Omegam]')
ylabel('Tiefe z [m]')
title('Darstellung mit der Tiefe')
hline1 = line(NaN,NaN, 'LineWidth', 2, 'LineStyle', '-', ...
    'Color', [0.85 0.85 0.85]);
hline2 = line(NaN,NaN, 'LineWidth', 2, 'LineStyle', '-', ...
    'Color', [0.4 0.4 0.4]);
hline3 = line(NaN,NaN, 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r');
legend([hline2, hline1, hline3], 'Quantil 25-75', 'Quantil 2.5-97.5', ...
    '\rho [\Omegam]')
grid on

%% Time display
et = toc;
etchar = sprintf('Elapsed time: %.2f sec. \n', et);
fprintf(etchar);

%% Save option

% save('variables.mat')
% rhoan = Het.rhoa_n;
% save('rhoa_n1D.mat','rhoan'); 

