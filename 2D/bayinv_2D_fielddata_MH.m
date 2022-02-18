% DC 2D problem for Checkerboard with MH algorithm
% Bayesian Inversion
% Author: Kristina Backes
% 2022

close all; clear; 
tic

%% Define problem params.

% Choose boundary condition type.
bc_type = 'dirichlet';

% Set up fwp.
% % ref_steps = 1;
FE_order = 1;
solver_type = util.pick(1, 'backslash', 'mumps', 'pcg_amg');

% Load measurement configuration.
% Poldipole
% [dc_info, d_obs] = app_dc.data.get_obs('name', 'erzp_2021.dat', ...
%      'verbosity', true);

% Dipoledipole
% [dc_info, d_obs] = app_dc.data.get_obs('name', 'erzd_2021.dat', ...
%    'verbosity', true);

% Wenner 
% [dc_info, d_obs] = app_dc.data.get_obs('name', 'erzw_2021.dat', ...
%     'verbosity', true);

% All data
[dc_info, d_obs] = app_dc.data.get_obs('name', 'all_2021.dat', ...
    'verbosity', true);

%% Set up mesh, bloc edge length 1m
% % %edge length
% % w = 1;
% % 
% % % distance ground surface to first bloc
% % di = 1; % 
% %
% % blockx1 = (0.5:w:19).';
% % blockx6 = 5;
% %
% % block0 = [-5, -1, 30, 0.99];
% % block1 = [blockx1, -di-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block2 = [blockx1, -di-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block3 = [blockx1, -di-3*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block4 = [blockx1, -di-4*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block5 = [blockx1, -di-5*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block6 = [blockx1, -di-6*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block7 = [blockx1, -di-7*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block8 = [blockx1, -di-8*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block9 = [blockx1, -di-9*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % 
% % block10 = [blockx6, -di-20*w+zeros(size(blockx6, 1), 1), 10+zeros(size(blockx6, 1), 1) ...
% %     5*w+zeros(size(blockx6, 1), 1) ];
% % 
% % block = [block0; block1; block2; block3; block4; block5; block6; block7; block8; ...
% %     block9; block10];

%% Set up mesh, bloc edge length 2m
% %edge length
% w = 2;
% 
% % distance ground surface to first bloc
% di = 1; % 
%
% blockx1 = (0:w:19).';
% blockx6 = 5;
% 
% block0 = [-5, -1, 30, 0.99];
% block1 = [blockx1, -di-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block2 = [blockx1, -di-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block3 = [blockx1, -di-3*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block4 = [blockx1, -di-4*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block9 = [blockx1, -di-10*1+zeros(size(blockx1, 1), 1), 2+zeros(size(blockx1, 1), 2)];
% block10 = [blockx6, -di-20*w+zeros(size(blockx6, 1), 1), 10+zeros(size(blockx6, 1), 1) ...
%     5*w+zeros(size(blockx6, 1), 1) ];
% 
% block = [block0; block1; block2; block3; block4; block9; block10];

%% Set up mesh, bloc edge length 0.5m
% %edge length
% w = 0.5;
% 
% % distance ground surface to first bloc
% di = 1; % 
%
% blockx1 = (0.5:w:19).';
% blockx6 = 5;
% 
% block0 = [-5, -1, 30, 0.99];
% block1 = [blockx1, -di-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block2 = [blockx1, -di-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block3 = [blockx1, -di-3*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block4 = [blockx1, -di-4*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block5 = [blockx1, -di-5*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block6 = [blockx1, -di-6*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block7 = [blockx1, -di-7*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block8 = [blockx1, -di-8*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block9 = [blockx1, -di-9*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block10 = [blockx1, -di-10*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block11 = [blockx1, -di-11*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block12 = [blockx1, -di-12*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block13 = [blockx1, -di-13*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block14 = [blockx1, -di-14*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block15 = [blockx1, -di-15*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block16 = [blockx1, -di-16*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block17 = [blockx1, -di-17*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block18 = [blockx1, -di-18*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block20 = [blockx6, -di-40*w+zeros(size(blockx6, 1), 1), 10+zeros(size(blockx6, 1), 1) ...
%     5*w+zeros(size(blockx6, 1), 1) ];
% 
% block = [block0; block1; block2; block3; block4; block5; block6; block7; block8; block9; block10; ...
%     block11; block12; block13; block14; block15; block16; block17; block18; block20];

%% Set up mesh, bloc edge length 1.5m
% %edge length
% w = 1.5;
% 
% % distance ground surface to first bloc
% di = 1; % 
%
% blockx1 = (0.25:w:19).';
% blockx6 = 5;
% 
% block0 = [-5, -1, 30, 0.99];
% block1 = [blockx1, -di-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block2 = [blockx1, -di-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block3 = [blockx1, -di-3*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block4 = [blockx1, -di-4*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block5 = [blockx1, -di-5*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block6 = [blockx1, -di-6*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block20 = [blockx6, -di-20*w+zeros(size(blockx6, 1), 1), 10+zeros(size(blockx6, 1), 1) ...
%     5*w+zeros(size(blockx6, 1), 1) ];
% 
% block = [block0; block1; block2; block3; block4; block5; block6; block20;];

%% Set up mesh, bloc edge length increase from 1m to 3m
% distance ground surface to first bloc
di = 1; % 

% edge length
w = 1;
blockx1 = (0.5:w:19).';
w3 = 2;
blockx2 = (1:w3:18).';
w4 = 3;
blockx4 = (1:w4:18).';

blockx6 = 5;

block0 = [-5, -1, 30, 0.99];
block1 = [blockx1, -di-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
block2 = [blockx1, -di-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
block4 = [blockx2, -di-4*w+zeros(size(blockx2, 1), 1), w3+zeros(size(blockx2, 1), 2)];
block6 = [blockx2, -di-6*w+zeros(size(blockx2, 1), 1), w3+zeros(size(blockx2, 1), 2)];
block9 = [blockx4, -di-9*w+zeros(size(blockx4, 1), 1), w4+zeros(size(blockx4, 1), 2)];

block10 = [blockx6, -di-20*w+zeros(size(blockx6, 1), 1), 10+zeros(size(blockx6, 1), 1) ...
    5*w+zeros(size(blockx6, 1), 1) ];

block = [block0; block1; block2; block4; block6; block9; block10];

%% Meshing
domain_r = 5e4;

[mesh, pn, pm, fm, cm] = meshing.generate_checkerboard2D(...
    'block', block, ...
    'point', dc_info.survey.ele, ...
    'keep_files', true, ...
    'domain_r', domain_r, ...
    'marker', [-1, 0, 1, 2]);

%% include anomaly
% Set up boundary conditions.
bc_facet_tag = pn{2}{1}(ismember(pn{2}{2}, {'subsurface'}));
bc_value = @(x) 0;
bnd_info = {fm, {bc_facet_tag}, {bc_value}};
bnd_summary = {bc_type, bnd_info};

%% Assembling and solving.

% Assemble fwd problem.
sol = app_dc.fwd.assemble('2.5D', mesh, FE_order, bnd_summary, dc_info);

% Assign measurements to noised data
rhoa_n = d_obs;

% Number of measurements
nn = length(rhoa_n);

%% Plot model and pseudosection

% number of blocs
len_split = length(block)+1; 

figure(41)
plot_model_check(mesh, pm, fm, cm, 'Checkerboard')
hold on
text(block(1,1)+5, block(1,2)+0.5, ""+(1+1), 'Fontsize', 10, 'color', [0 0 0])
for k = 2:39
    text(block(k,1), block(k,2)+0.5, ""+(k+1), 'Fontsize', 10, 'color', [0 0 0])
end
for k = 40:len_split-2
    text(block(k,1), block(k,2)+0.5, ""+(k+1), 'Fontsize', 14, 'color', [0 0 0])
end
text(block(len_split-1,1), block(len_split-1,2)+2.5, ""+(len_split-1+1), 'Fontsize', 14, 'color', [0 0 0])
hold off
set(gca, 'Fontsize', 16)
setcolorbartitle('Box number')
xlim([-2 22])
ylim([-20 2])

figure(42)
app_dc.data.plot_pseudosection(dc_info.survey, rhoa_n)
title('Pseudosektion - Pol-Dipol')
set(gca, 'Fontsize', 16)
setcolorbartitle('[\Omegam]')

figure(43)
plot_model_check(mesh, pm, fm, cm, 'Checkerboard')
set(gca, 'Fontsize', 16)
setcolorbartitle('Box number')
xlim([-2 22])
ylim([-20 2])

%% Define MCMC

% number iterations
n = 39000; 

% burnin amount
burn = 10000;

% standard deviation of normal proposal density or step size
% of conductivity value
s = 0.05;

% Variance
sigma = 0.1; 

% Auxiliary variable 
mittelA = 0;

% Random walk of resistivity
target_leitf = zeros(len_split, n); 
target_leitf(:,1) = 1/500 * ones(len_split, 1);

% Declare variables
proposed_log_rho = zeros(len_split, 1);
result.layer = zeros(n-1, len_split);
result.stat = zeros(n-1, 2);

%% Startmodell plot
start_rhoa = log(1./target_leitf(:,1));
% Fixed background resistivity
start_rhoa(1) = log(500);
start_rhoa(end) = log(500);

% First layer
start_rhoa(2) = log(500); % log(200)

param_start = (cm == unique(cm).') * start_rhoa;

figure(52)
plot_model(mesh, 1./exp(param_start), 'Model Prior', dc_info)
set(gca, 'Fontsize', 16)
setcolorbar('[\Omegam]')
xlim([-10 30])
ylim([-15 2])

%% Metropolis-Hastings-Algorithm

% Set acceptance to true
hlp_time = 1;

% MCMC loop
for i = 2:n
    
    % Assignement of previous position in MCMC chain
    current_log_rho = log(1./target_leitf(:,i-1));
    
    % Generating new candidate for random walk
    for k = 1:len_split 
        proposed_log_rho(k) = current_log_rho(k) + s*randn(1,1);
    end
    
    % Fixed background resistivity
    proposed_log_rho(1) = log(500);
    proposed_log_rho(end) = log(500);
    
    % Fixed first layer
    %proposed_log_rho(2) = log(100); % log(200)
    
    % Set up parameter space candidate
    param_prop = (cm == unique(cm).') * proposed_log_rho;
    
    % Forward operation candidate
    fwd_propr = app_dc.fwd.solve(sol, 1./exp(param_prop), ...
        solver_type, dc_info); 
    
    % Prior candidate
    prior_prop = normpdf(proposed_log_rho, 4.8, 0.9);
    prior_prop(2) = normpdf(proposed_log_rho(2), 4.8, 0.2);
         
    % Likelihood candidate
    % Relative Log-Like
    log_like_prop = - sum(log(fwd_propr)) - (1/(2*sigma^2)) * ...
        sum((rhoa_n - fwd_propr).^2 ./ fwd_propr.^2);
    
    % If previous step was accepted
    if hlp_time == 1
        
        % Set up parameter space previous value i-1
        param_curr = (cm == unique(cm).') * current_log_rho;
    
        % Forward operation previous value i-1
        fwd_currr = app_dc.fwd.solve(sol, 1./exp(param_curr), ...
            solver_type, dc_info);
        
        % Prior previous value i-1
        prior_curr = normpdf(current_log_rho, 4.8, 0.9);
        prior_curr(2) = normpdf(current_log_rho(2), 4.8, 0.2);
        
        % Likelihood previous value i-1
        % Relative Log-Like
        log_like_curr = - sum(log(fwd_currr)) - (1/(2*sigma^2)) * ...
            sum((rhoa_n - fwd_currr).^2 ./ fwd_currr.^2); 
    end
   
    % Posterior density
    post_dens = sum(log(prior_prop(:))) + log_like_prop;
        
    % Probability of move
    A = exp(post_dens - sum(log(prior_curr)) - log_like_curr);
    
    % Generate uniformly distributed random variable
    u = rand(1,1);
    
    % Acceptance query
    if u <= min(A, 1)
        % Random walk
        target_leitf(:,i) = 1./exp(proposed_log_rho);
        
        % Results
        result.layer(i-1,:) = exp(proposed_log_rho);
        result.stat(i-1,1) = log_like_prop;
        result.stat(i-1,2) = post_dens;
        
        % Auxiliary variable
        hlp_time = 1;
        mittelA = mittelA + 1;
        
    else
        % Random walk
        target_leitf(:,i) = 1./exp(current_log_rho); 
        
        % Auxiliary variable
        hlp_time = 0;
    end
    
    % Progress
    if i == round(0.25*n)
        fprintf('Fortschritt: 25 %% \n');
    end
    if i == round(0.5*n)
        fprintf('Fortschritt: 50 %% \n');
    end
    if i == round(0.75*n)
        fprintf('Fortschritt: 75 %% \n');
    end
end

% Acceptance probability of random walk
akwrs = (mittelA/n) * 100; 

%% Estimator 1 - MAP

% Sort for maximal posteriori density/likelihood
[~, s_idx] = sort(result.stat(:, 1), 'descend'); 
result.stat = result.stat(s_idx, :); 
result.layer = result.layer(s_idx, :); 

idx = find(result.stat(:,1), 1); 
result.stat = result.stat(idx:end, :); 
result.layer = result.layer(idx:end, :);

% Computation of the map_n most probable modell answers
map_n = 20;
model_map = zeros(len_split, map_n);
for l = 1:map_n
    model_map(:,l) = result.layer(l,:);
end

% Statistics
model_mean = mean(model_map,2);
model_var = var(model_map,0,2);
model_std = std(model_map,0,2);
model_stat = [model_mean model_std];

% Quantile 
rho_quant = quantile(model_map,[0.025 0.25 0.50 0.75 0.975],2);

%% Estimator 2 - Mean of Markov Chain
% Statistics
model_mean2 = mean(1./target_leitf(:,burn:end),2);
model_var2 = var(1./target_leitf(:,burn:end),0,2);
model_std2 = std(1./target_leitf(:,burn:end),0,2);
model_stat2 = [model_mean2 model_std2];

% Quantile 
rho_quant2 = quantile(1./target_leitf(:,burn:end), ...
    [0.025 0.25 0.50 0.75 0.975],2);

model_std_prozentual = model_std2 ./ model_mean2;

%% Plot Random Walk

pax = min(len_split, 4*6);
maxxi = max(1./target_leitf, [], 'all');

figure(11)
for m = 1:pax
    subplot(4,6,m)
    hold on
    plot(1:n, 1./target_leitf(m+1,:), 'ko'); 
    plot([burn, burn], [1, maxxi], 'r-', 'linewidth', 1.5);
    set(gca, 'YScale', 'log')
    title("Block " + (m+1)); 
    xlabel('Index')
    ylabel('Zustand \rho [\Omegam]')
    xlim([0, n])
    ylim([min(1./target_leitf(m+1,:)), max(1./target_leitf(m+1,:))])
    grid on
end

%% Priori-distribution 

for i = 1:n
    prio_hist(i) = lognrnd(4.8, 0.9); 
end

%% Comparison Priori & Posteriori

for k = 1:len_split
    [f(k,:), xi(k,:)] = ksdensity(1./target_leitf(k,:), (1:2:4000));
end

prior_histks = ksdensity(prio_hist, (1:2:4000));

figure(102)
for i = 1:(2*5)
    subplot(2,5,i)
    set(gca, 'Fontsize', 12)
    hold on
    plot(xi(1,:), prior_histks, 'LineWidth', 2.5)
    plot(xi(1,:), f(i+50,:), 'LineWidth', 2.5)
    set(gca, 'XScale', 'log')
    title("Block" + (i+50))
    xlim([10 4000])
    xticks(logspace(1,4,4))
    xlabel('\rho [\Omegam]')
    if i == 1 || i == 6
        ylabel("PDF")
    end
    if i == 1
        legend("\mu_0", "\pi", 'location', 'northeast' ,'Fontsize',13)
    end
    grid on
    ax = gca;
    ax.YAxis.Exponent = 0;
end


%% Error plots Mean of Chain

% Set up parameter space with model mean
param_result = (cm == unique(cm).') * model_stat2(:,1);

% Set up parameter space with model mean
param_result_std = (cm == unique(cm).') * model_stat2(:,2);

% Forward operation
fwd_result = app_dc.fwd.solve(sol, 1./param_result, solver_type, dc_info);

% Error
err_rhoa = abs(fwd_result - rhoa_n);

% Colorbar boundary
comin = min(param_result, [], 'all');
comax = max(param_result, [], 'all');

% Plots seudosection
figure(21)
app_dc.data.plot_pseudosection(dc_info.survey, fwd_result)
title('Result Mean')
set(gca,'ColorScale','log')
set(gca, 'Fontsize', 16)
setcolorbartitle(' \rho_a [\Omegam]')

figure(22)
app_dc.data.plot_pseudosection(dc_info.survey, err_rhoa)
title('Absoluter Fehler')
set(gca, 'Fontsize', 16)
setcolorbartitle('\rho_a [\Omegam]')

%% Plot model outcome Mean
figure(23)
plot_model(mesh, 1./param_result, 'Model Mean', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
caxis([comin comax])
setcolorbar('\rho [\Omegam]', comin, comax)
xlim([-3 23])
ylim([-12 2])

%%
figure(24)
plot_model(mesh, 1./param_result_std, 'Model Mean Std', dc_info)
set(gca, 'Fontsize', 16)
%set(gca,'ColorScale','log')
setcolorbartitle('\rho [\Omegam]')
xlim([-3 23])
ylim([-12 2])

%% Data misfit

figure(101)
plot_data_misfit(fwd_result, d_obs)
set(gca, 'Fontsize', 16)
grid on

%% Error plots MAP

% Set up parameter space
param_result_map = (cm == unique(cm).') * model_stat(:,1);

% Set up parameter space with model mean
param_result_map_std = (cm == unique(cm).') * model_stat(:,2);

% Forward operation
fwd_result_map = app_dc.fwd.solve(sol, 1./param_result_map, solver_type, dc_info);

% Error
err_rhoa_map = abs(fwd_result_map - rhoa_n);

% Plots pseudosection
figure(31)
set(gca, 'Fontsize', 13)
app_dc.data.plot_pseudosection(dc_info.survey, fwd_result_map)
title('Result MCMC 1 MAP')

figure(32)
set(gca, 'Fontsize', 13)
app_dc.data.plot_pseudosection(dc_info.survey, err_rhoa_map)
title('Error Rhoa 1 MAP')

%% Plot model outcome MAP
figure(33)
plot_model(mesh, 1./param_result_map, 'Model MAP', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
caxis([comin comax])
xlim([-3 23])
ylim([-12 2])
setcolorbar('\rho [\Omegam]', comin, comax)

figure(34)
plot_model(mesh, 1./param_result_map_std, 'Model MAP Std', dc_info)
set(gca, 'Fontsize', 16)
setcolorbartitle('\rho [\Omegam]')
xlim([-3 23])
ylim([-12 2])

%% Resultate für MAP, most probable

% Set up parameter space
param_result_map1 = (cm == unique(cm).') * result.layer(1,:)';

% Forward operation
fwd_result_map1 = app_dc.fwd.solve(sol, 1./param_result_map1, solver_type, dc_info);

% Error
err_rhoa_map1 = abs(fwd_result_map1 - rhoa_n);

% Plots
figure(61)
plot_model(mesh, 1./param_result_map1 , 'MAP 1 Modell', dc_info)
set(gca, 'Fontsize', 16)
xlim([-5 25])
ylim([-15 2])

%% Time display
et = toc;
etchar = sprintf('Elapsed time: %.2f sec. \n', et);
fprintf(etchar);

%% Save option

% save('15_mh_variables_alldata.mat');

%% Helper functions

function plot_model(mesh, sigma, str, dc_info)

    meshing.plot_mesh(mesh, 'cell_markers', 1./sigma);
    set(gca, 'Fontsize', 13)
    ylabel('z [m]');
    title(str);
    xlabel('x [m]');
    hold on
        plot(dc_info.survey.ele(:, 1), dc_info.survey.ele(:, 2), 'c.');
%         meshing.plot_mesh(mesh)
    hold off
    drawnow;
end

function plot_model_check(mesh, pm, fm, cm, str)

    meshing.plot_mesh(mesh, 'vertex_markers', pm, ...
                        'facet_markers', fm, ...
                        'cell_markers', cm);
    set(gca, 'Fontsize', 13)
    ylabel('z [m]');
    title(str);
    xlabel('x [m]');
%     hold on
%         meshing.plot_mesh(mesh)
%     hold off
    drawnow;
end

function plot_data_misfit(d_fwd, d_obs)

    n_data = length(d_obs);
    figure(101)
        plot(1:n_data, d_fwd, '.r', ...
             1:n_data, d_obs, 'ob');
        set(gca, 'YScale', 'log')
        legend('d_{i}', 'd_{obs}');
        title('Messwertevergleich');
        xlabel('Nummer der Messung');
        ylabel('\rho_a [\Omegam]');
end

function setcolorbartitle(str)
    colorcet('D1')
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = str;
    set(colorTitleHandle ,'String',titleString);
    hcb.Ruler.Exponent = 0;
end

function setcolorbar(str, comin, comax)
    colorcet('D1')
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = str;
    set(colorTitleHandle ,'String',titleString);
    set(hcb,'YTick', round(logspace(log10(comin), log10(comax),7),0))
    hcb.Ruler.Exponent = 0;
end