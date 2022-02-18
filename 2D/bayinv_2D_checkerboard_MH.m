% DC 2D problem for Checkerboard
% Bayesian Inversion
% Author: Kristina Backes
% 2022

close all; clear; %clc;
tic

%% Define problem params.

% Analysis slim or thorough
interim_result = util.pick(1, true, false);

% Choose boundary condition type.
bc_type = 'dirichlet';

% Set up fwp.
ref_steps = 0;
FE_order = 1;
solver_type = util.pick(2, 'backslash', 'mumps', 'pcg_amg');

% Set up domain.
% Set electrode positions at top of halfspace.
topo_pos = [(-100:100).', 0*(-100:100).'];

% Define conductivity of halfspace/background
rho.back = 1000;
rho.ano = 400;
domain_val = 1 / rho.back; 

% Define electrode positions.
ixd_ele_in_topo = topo_pos(:,1) > -25 & topo_pos(:,1) < 25;
ele_pos = topo_pos(ixd_ele_in_topo, :);

% Create a measurement configuration.                           
dc_info = app_dc.data.create_survey(ele_pos, 'poledipole');
dc_info.survey.data_type = 'rhoa';
 
%% Set up mesh, edge length 3m
% w = 3;
% blockx1 = (-14:w:12).';
% 
% block1 = [blockx1, -1-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block2 = [blockx1, -1-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block3 = [blockx1, -1-3*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block4 = [blockx1, -1-4*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block6 = [0, -30, 5, 5];
% block = [block1; block2; block3; block4; block6];

%% Set up mesh, edge length 6m
% w = 6;
% blockx1 = (-14:w:12).';
% 
% block1 = [blockx1, -1-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block2 = [blockx1, -1-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% block6 = [0, -30, 5, 5];
% block = [block1; block2; block6];

%% Set up mesh, edge length 1m
% % w = 1;
% % blockx1 = (-14:w:12).';
% % block = [blockx1, -1-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% % block6 = [0, -30, 5, 5];
% % 
% % di = 1;
% % 
% % for k = 1:8
% %     B = [blockx1, -1-di-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
% %     di = di+1;
% %     block = vertcat(block, B);
% % end
% % 
% % block = vertcat(block, block6);

%% Set up parameter vector

% domain_val_vec = (1/rho.ano) * ones(length(block) + 1, 1); % 1/100 -> Anomaly 100 Ohmm
% domain_val_vec(1:2:end) = domain_val; % Hintergrund 800 Ohmm
% domain_val_vec(end) = domain_val; % letzter Block
% 
% param_fine = (cm_fine == unique(cm_fine).') * domain_val_vec;

%%  Set up mesh, edge length 3m, anomaly 6m

w = 3;
blockx1 = (-14:w:14).';

block1 = [blockx1, -1-w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
block2 = [blockx1, -1-2*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
block3 = [blockx1, -1-3*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
block4 = [blockx1, -1-4*w+zeros(size(blockx1, 1), 1), w+zeros(size(blockx1, 1), 2)];
block6 = [0, -30, 5, 5];
block = [block1; block2; block3; block4; block6];

domain_r = 5e4;

[mesh_fine, pn_fine, pm_fine, fm_fine, cm_fine] = meshing.generate_checkerboard2D(...
    'block', block, ...
    'point', dc_info.survey.ele, ...
    'keep_files', true, ...
    'ref', ref_steps, ...
    'domain_r', domain_r, ...
    'marker', [-1, 0, 1, 2]);

figure(20)
plot_model_check(mesh_fine, pm_fine, fm_fine, cm_fine, 'Checkerboard')

domain_val_vec = domain_val * ones(4, length(blockx1)); 

for k = 1:2
    domain_val_vec(k,1:4:end) = 1/rho.ano; 
    domain_val_vec(k,2:4:end) = 1/rho.ano; 
end
for k = 3:4
    domain_val_vec(k,3:4:end) = 1/rho.ano;
    domain_val_vec(k,4:4:end) = 1/rho.ano; 
end

domain_val_vec = reshape(domain_val_vec', [], 1);

domain_val_vec = vertcat(domain_val_vec, domain_val); % letzter Block
domain_val_vec = vertcat(domain_val, domain_val_vec); % erster Block

param_fine = (cm_fine == unique(cm_fine).') * domain_val_vec;

%% Set up boundary conditions.
bc_facet_tag = pn_fine{2}{1}(ismember(pn_fine{2}{2}, {'subsurface'}));
bc_value = @(x) 0;
bnd_info = {fm_fine, {bc_facet_tag}, {bc_value}};
bnd_summary = {bc_type, bnd_info};

%% Assembling and solving.

% Assemble fwd problem.
sol_fine = app_dc.fwd.assemble('2.5D', mesh_fine, FE_order, bnd_summary, dc_info);

% Solve fwd problem.
rhoa = app_dc.fwd.solve(sol_fine, param_fine, solver_type, dc_info);

% noise type: 
% 1: additive noise, 2: relative noise, 3: superposition 
noise_type = 2;

switch noise_type
    case 1
        % Additive noise
        noise = 3;
        rhoa_n = rhoa + noise .* randn(size(rhoa));
        sigma = noise; % Standardabweichung
    case 2
        % Relative noise
        noise = 0.03;
        rhoa_n = rhoa .* (1 + noise * randn(size(rhoa)) );
        sigma = noise; % Varianz
    case 3
        % Superposition relativ and additive
        noise = 0.03;
        noise_g = 0.003;
        noise_geo = noise_g * AB + zeros(size(AB));
        rhoa_nn = rhoa .* (1 + noise * randn(size(rhoa)) );
        rhoa_n = rhoa_nn + noise_geo' .* randn(size(rhoa));
        sigma = noise; % Varianz
end

rhoa_n_log = log(rhoa_n);

% number of observations
nn = length(rhoa);

%% Interim result plots

if interim_result
    figure(41)
    set(gca, 'Fontsize', 12)
    app_dc.data.plot_pseudosection(dc_info.survey, rhoa)
    title('Resistivity [\Omegam]')

    figure(42)
    set(gca, 'Fontsize', 12)
    plot_model(mesh_fine, param_fine, 'Model with anomaly', dc_info)
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = '[\Omegam]';
    set(colorTitleHandle ,'String',titleString);
    xlim([-20 20])
    ylim([-30 5])
end

%% Load data

% All data
% % [dc_info, rhoa_n] = app_dc.data.get_obs('name', '/home/kris/MA/WichtigeDatensaetze/Checkerboard_BI_backes.dat', ...
% %     'verbosity', true);
% % 
% % rhoa_n_log = log(rhoa_n);
% % nn = length(rhoa_n);
% % noise = 0.03;
% % sigma = noise; % Varianz

% larger mesh for inversion
[mesh, pn, pm, fm, cm] = meshing.generate_checkerboard2D(...
    'block', block, ...
    'point', dc_info.survey.ele, ...
    'keep_files', true, ...
    'domain_r', domain_r, ...
    'marker', [-1, 0, 1, 2]);

% Set up boundary conditions.
bc_facet_tag = pn{2}{1}(ismember(pn{2}{2}, {'subsurface'}));
bc_value = @(x) 0;
bnd_info = {fm, {bc_facet_tag}, {bc_value}};
bnd_summary = {bc_type, bnd_info};

% Assemble fwd problem.
sol = app_dc.fwd.assemble('2.5D', mesh, FE_order, bnd_summary, dc_info);

%% Define MCMC

% number iterations
n = 35000; 

% burnin amount
burn = 10000;

% standard deviation of normal proposal density or step size
s = 0.03;

% Auxiliary variable
len_split = length(block)+1; 
mittelA = 0;

% Random walk of resistivity
target_leitf = zeros(len_split, n);
target_leitf(:,1) = 1/rho.back * ones(len_split, 1);

% Declare variables
proposed_log_rho = zeros(len_split, 1);
result.layer = zeros(n-1, len_split);
result.stat = zeros(n-1, 2);

%% MH-Algorithm

% Set acceptance to true
hlp_time = 1;

% MCMC Loop
for i = 2:n
    
    % Assignement of previous position in MCMC chain
    current_log_rho = log(1./target_leitf(:,i-1));
    
    % Generating new candidate for random walk
    for k = 1:len_split % Jedes rho unterschiedlich stören
        proposed_log_rho(k) = current_log_rho(k) + s*randn(1,1);
    end
    
    % Fixed surroundings
%     proposed_log_rho(1) = log(rho.back);

    % Connecting surroundings with numeric extra bloc
    proposed_log_rho(end) = proposed_log_rho(1);
    
    % Set up parameter space candidate
    param_prop = (cm == unique(cm).') * proposed_log_rho;
    
    % Forward operation candidate
    fwd_propr = app_dc.fwd.solve(sol, 1./exp(param_prop), ...
        solver_type, dc_info);
    fwd_propr_log = log(fwd_propr);
    
    % Prior candidate
     prior_prop = normpdf(proposed_log_rho, 4.2, 0.8);
%     prior_prop = normpdf(proposed_log_rho, 5.8, 0.5);
         
    % Likelihood candidate
    switch noise_type
        case 1
            % Additive Log-Like
            log_like_prop = - (1/(2*sigma^2)) * ...
                sum((rhoa_n_log - fwd_propr_log).^2);
        case 2
            % Relative Log-Like
            log_like_prop = - sum(log(fwd_propr_log)) - (1/(2*sigma^2)) * ...
                sum((rhoa_n_log - fwd_propr_log).^2 ./ fwd_propr_log.^2);
        case 3
            % Superposition
            log_like_prop = - sum(log(hlp_prop * noise + noise_geo')) - ...
                0.5 * sum((rhoa_n_log - fwd_propr_log).^2 ./ ...
                (fwd_propr_log.^2 * noise^2 + noise_geo'.^2));
    end
    
    % If previous step was accepted
    if hlp_time == 1
        
        % Set up parameter space previous value i-1
        param_curr = (cm == unique(cm).') * current_log_rho;
        
        % Forward operation previous value i-1
        fwd_currr = app_dc.fwd.solve(sol, 1./exp(param_curr), ...
            solver_type, dc_info); % rhoa für current_x
        fwd_currr_log = log(fwd_currr);
        
        % Likelihood previous value i-1
        switch noise_type
            case 1
                % Additive Log-Like
                log_like_curr = - (1/(2*sigma^2)) * ...
                    sum((rrhoa_n_log - fwd_currr_log).^2);
            case 2
                % Relative Log-Like
                log_like_curr = - sum(log(fwd_currr_log)) - (1/(2*sigma^2)) * ...
                    sum((rhoa_n_log - fwd_currr_log).^2 ./ fwd_currr_log.^2); 
            case 3
                % Superposition
                log_like_curr = - sum(log(fwd_currr_log * noise + noise_geo')) - ...
                    0.5 * sum((rhoa_n_log - fwd_currr_log).^2 ./ ...
                    (ffwd_currr_log.^2 * noise^2 + noise_geo'.^2));
        end
        
        % Prior previous value i-1
         prior_curr = normpdf(current_log_rho, 4.2, 0.8);
%        prior_curr = normpdf(current_log_rho, 5.8, 0.5);
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
        target_leitf(:,i) = 1./exp(proposed_log_rho); %accept proposal as new
        
        % Results
        result.layer(i-1,:) = exp(proposed_log_rho);
        result.stat(i-1,1) = log_like_prop;
        result.stat(i-1,2) = post_dens ;
        
        % Auxiliary variable
        hlp_time = 1;
        mittelA = mittelA + 1; % akzeptiert neuen Move
        
    else
        % Random walk
        target_leitf(:,i) = 1./exp(current_log_rho); % and set old to be the new
        
        % Auxiliary variable
        hlp_time = 0;
    end
    
    % Progress
    if i == round(0.25*n)
        fprintf('Fortschritt: 25 %%\n');
    end
    if i == round(0.5*n)
        fprintf('Fortschritt: 50 %%\n');
    end
    if i == round(0.75*n)
        fprintf('Fortschritt: 75 %%\n');
    end
end

% Acceptance probability of random walk
akwrs = (mittelA/n) * 100; 

%% Estimator 1 - MAP

% Sort for maximal posteriori density = maximal likeliihood
[~, s_idx] = sort(result.stat(:, end), 'descend'); 
result.stat = result.stat(s_idx, :); 
result.layer = result.layer(s_idx, :); 

idx = find(result.stat(:,1), 1); 
result.stat = result.stat(idx:end, :); 
result.layer = result.layer(idx:end, :);

% Computation of the map_n most probable modell answers
map_n = 20;
model_map = zeros(len_split, map_n);
for i = 1:map_n
    model_map(:,i) = result.layer(i,:);
end

% Statistics
model_mean = mean(model_map,2);
model_var = var(model_map,0,2);
model_std = std(model_map,0,2);
model_stat = [1./model_mean 1./model_std];

% Quantile 
rho_quant = quantile(model_map,[0.025 0.25 0.50 0.75 0.975],2);

%% Estimator 2 - Mean of Markov Chain
% Statistics
model_mean2 = mean(1./target_leitf(:,burn:end),2);
model_var2 = var(1./target_leitf(:,burn:end),0,2);
model_std2 = std(1./target_leitf(:,burn:end),0,2);
model_stat2 = [model_mean2 model_std2];

% Quantile 
rho_quant2 = quantile(1./target_leitf(:,burn:end),...
    [0.025 0.25 0.50 0.75 0.975],2);

%% Plot Random Walk

figure(1)
min_len = min(len_split, 30);
for m = 1:min_len
    subplot(5,6,m)
    hold on
    plot(1:n, 1./target_leitf(m,:), 'ko'); 
    plot(1:n, ones(1,n)*(1./domain_val_vec(m)), 'LineWidth', 2);
    title("MCMC-Walk layer " + m + " (" + 1./domain_val_vec(m) + "\Omegam)")
    xlabel('Index')
    ylabel('Besuchte Orte')
    xlim([0,n])
end

%% Error plots

param = (cm == unique(cm).') * domain_val_vec;

% Set up parameter space
param_result = (cm == unique(cm).') * model_stat2(:,1);
param_result_std = (cm == unique(cm).') * model_stat2(:,2);

% Forward operation
fwd_result = app_dc.fwd.solve(sol, 1./param_result, solver_type, dc_info);

% Error
err_rhoa = abs(fwd_result - rhoa_n);

cominr = min([rhoa_n; fwd_result], [], 'all');
comaxr = max([rhoa_n; fwd_result], [], 'all');

% Plots
figure(21)
set(gca, 'Fontsize', 15)
app_dc.data.plot_pseudosection(dc_info.survey, rhoa_n)
title('Messwert \rho_a')
set(gca,'ColorScale','log')
setcolorbartitle('\rho_a [\Omegam]')

figure(22)
set(gca, 'Fontsize', 15)
app_dc.data.plot_pseudosection(dc_info.survey, fwd_result)
title('Result Mean')
set(gca,'ColorScale','log')
setcolorbartitle('\rho_a [\Omegam]')

figure(23)
set(gca, 'Fontsize', 15)
app_dc.data.plot_pseudosection(dc_info.survey, err_rhoa)
title('Absoluter Fehler')
setcolorbartitle('\rho_a [\Omegam]')

%% Plot model with elektrode
paramr = 1./param;

comin = min([paramr; param_result], [], 'all');
comax = max([paramr; param_result], [], 'all');

figure(31)
plot_model(mesh, 1./paramr, 'Model', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
setcolorbar('\rho [\Omegam]', comin, comax)
caxis([300 1000])
xlim([-20 20])
ylim([-15 2])

figure(32)
plot_model(mesh, 1./param_result, 'Result Mean',dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
setcolorbar('\rho [\Omegam]', comin, comax)
caxis([300 1000])
xlim([-20 20])
ylim([-15 2])

figure(33)
plot_model(mesh, 1./param_result_std, 'Model Mean Std', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
setcolorbar('\rho [\Omegam]', min(param_result_std), max(param_result_std));
caxis([min(param_result_std) max(param_result_std)])
xlim([-20 20])
ylim([-15 2])

%% Time display
et = toc;
etchar = sprintf('Elapsed time: %.2f sec. \n', et);
fprintf(etchar);

%% Save option

% save('15_mh_variables_checker.mat')

% app_dc.data.export2BERT(dc_info.survey, 'Checkerboard_BI_pd2')

%% Help functions

function plot_model(mesh, sigma, str, dc_info)

    meshing.plot_mesh(mesh, 'cell_markers', 1./sigma);
    ylabel('z [m]');
    title(str);
    xlabel('x [m]');
    hold on
        plot(dc_info.survey.ele(:, 1), dc_info.survey.ele(:, 2), 'ko');
        meshing.plot_mesh(mesh)
    hold off
    drawnow;
end

function plot_model_check(mesh, pm, fm, cm, str)

    meshing.plot_mesh(mesh, 'vertex_markers', pm, ...
                        'facet_markers', fm, ...
                        'cell_markers', cm);
    ylabel('z [m]');
    title(str);
    xlabel('x [m]');
    hold on
        meshing.plot_mesh(mesh)
    hold off
    drawnow;
end

function setcolorbartitle(str)
    colorcet('D1')
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = str;
    set(colorTitleHandle ,'String',titleString);
end

function setcolorbar(str, comin, comax)
    colorcet('D1')
    hcb = colorbar;
    colorTitleHandle = get(hcb,'Title');
    titleString = str;
    set(colorTitleHandle ,'String',titleString);
    set(hcb,'YTick', round(logspace(log10(comin), log10(comax),6),0))
end