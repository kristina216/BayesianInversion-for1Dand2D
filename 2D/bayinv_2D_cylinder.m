% DC 2D problem for cylindric anomaly
% Bayesian Inversion
% Author: Kristina Backes
% 2022

close all; clear; %clc;
tic

%% Define problem params.

% Choose boundary condition type.
bc_type = 'dirichlet';

% Set up fwp.
ref_steps = 1;
FE_order = 1;
solver_type = util.pick(2, 'backslash', 'mumps', 'pcg_amg');

% Set up domain.
% Set electrode positions at top of halfspace.
topo_pos = [(-100:100).', 0*(-100:100).'];

domain_val = 1 / 800; 

% Define electrode positions.
ixd_ele_in_topo = topo_pos(:,1) > -25 & topo_pos(:,1) < 25;
ele_pos = topo_pos(ixd_ele_in_topo, :);

% Create a measurement configuration.
dc_info = app_dc.data.create_survey(ele_pos, 'poledipole');
dc_info.survey.data_type = 'rhoa';

% Set up mesh.
[mesh, cm, fm, pn] = meshing.generate_mesh2D('point', dc_info.survey.ele, ...
                                             'topo', ele_pos, ...
                                             'ref', ref_steps, ...
                                             'size_at_pt', 2, ...
                                             'domain_r', 1e4, ...
                                             'marker', [2, 1, -1]);

% Set up parameter vector.
param = (cm == unique(cm).') * domain_val;

%% include anomaly
midpoints = mesh.get_cell_centroids;

kreis = anomaly_triangle(4, -10, 8, midpoints); % x, y, r, midpoints
% Conductivity of anomly
param(kreis) = 1 / 300;

% Set up boundary conditions.
bc_facet_tag = pn{2}{1}(ismember(pn{2}{2}, {'bnd_earth'}));
bc_value = @(x) 0;
bnd_info = {fm, {bc_facet_tag}, {bc_value}};
bnd_summary = {bc_type, bnd_info};

%% Assembling and solving.

% Assemble fwd problem.
sol = app_dc.fwd.assemble('2.5D', mesh, FE_order, bnd_summary, ...
                          dc_info);

% Solve fwd problem.
rhoa = app_dc.fwd.solve(sol, param, solver_type, dc_info);

% noise type: 
% 1: additive noise, 2: relative noise, 3: superposition 
noise_type = 2;

switch noise_type
    case 1
        % Additive noise
        noise = 3;
        rhoa_n = rhoa + noise .* randn(size(rhoa));
        sigma = noise; 
    case 2
        % Relative noise
        noise = 0.03;
        rhoa_n = rhoa .* (1 + noise * randn(size(rhoa)) );
        sigma = noise; 
    case 3
        % Superposition relativ and additive
        noise = 0.03;
        noise_g = 0.003;
        noise_geo = noise_g * AB + zeros(size(AB));
        rhoa_nn = rhoa .* (1 + noise * randn(size(rhoa)) );
        rhoa_n = rhoa_nn + noise_geo' .* randn(size(rhoa));
        sigma = noise; 
end

rhoa_n_log = log(rhoa_n);

% Number of measurements
nn = length(rhoa);

%% Interim results

figure(1)
app_dc.data.plot_pseudosection(dc_info.survey, rhoa)
title('Rhoa')

figure(2)
plot_model(mesh, param, 'Halbraum', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
setcolorbartitle('\rho [\Omegam]')

%% Define MCMC

% number iterations
n = 35000; 

% burnin amount
burn = 10000;

% standard deviation of normal proposal density or step size
s = 0.05; 
s_geo = 0.05;

% Auxiliary variable 
len_split = 2;
len_geo = 3; 
mittelA = 0;

% Ranom walk resistivity
target_leitf = zeros(len_split, n); 
% Start model = [conductivity of half space, conductivity of anomaly];
target_leitf(:,1) = [1/500; 1/500]; 

% Random walk geometry
target_geo = zeros(len_geo,n);
% Start model = [x, z, r];
target_geo(:,1) = [-2, -13, 5]; % vorher -10, -10, 5

% Declare variables
proposed_log_rho = zeros(1, len_split);
proposed_geo = zeros(1, len_geo);
param_prop = ones(size(param));
param_curr = ones(size(param));

%% Probability Density Function

% Set acceptance to true
hlp_time = 1;

% Define PDFs for geometry factors
pdx = makedist('Uniform','lower',-18,'upper',18);
pdy = makedist('Uniform','lower',-18,'upper',-12);
pdr = makedist('Uniform','lower',4,'upper',10);

for i = 2:n
    
    % Assignement of previous position in MCMC chain
    current_log_rho = log(1./target_leitf(:,i-1));
    current_geo = target_geo(:,i-1);
    
    % Generating new candidate for random walk
    for k = 1:len_split
        proposed_log_rho(k) = current_log_rho(k) + s*randn(1,1);
    end
    for l = 1:len_geo
        proposed_geo(l) = current_geo(l) + s_geo*randn(1,1);
    end
    
    % Set up parameter space candidate
    param_prop(:) = proposed_log_rho(1);
    kreis_prop = anomaly_triangle(proposed_geo(1), proposed_geo(2), ...
        proposed_geo(3), midpoints);
    param_prop(kreis_prop) = proposed_log_rho(2);
    
    % Forward operation candidate
    fwd_propr = app_dc.fwd.solve(sol, 1./exp(param_prop), solver_type, dc_info); % rhoa für proposed_x
    fwd_propr_log = log(fwd_propr);

    % Prior candidate
    prior_prop = normpdf(proposed_log_rho, 5.5, 0.9);
    prior_prop(len_split+1) = pdf(pdx, proposed_geo(1)); 
    prior_prop(len_split+2) = pdf(pdy, proposed_geo(2)); 
    prior_prop(len_split+3) = pdf(pdr, proposed_geo(3)); 
    
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
            log_like_prop = - sum(log(fwd_propr_log * noise + noise_geo')) - ...
                0.5 * sum((rhoa_n_log  - fwd_propr_log).^2 ./ ...
                (fwd_propr_log.^2 * noise^2 + noise_geo'.^2));
    end
    
    % If previous step was accepted
    if hlp_time == 1  
        % Prior previous value i-1
        prior_curr = normpdf(current_log_rho, 5.5, 0.9);
        prior_curr(len_split+1) = pdf(pdx, current_geo(1)); 
        prior_curr(len_split+2) = pdf(pdy, current_geo(2)); 
        prior_curr(len_split+3) = pdf(pdr, current_geo(3));
        
        % Set up parameter previous value i-1
        param_curr(:) = current_log_rho(1);
        kreis_curr = anomaly_triangle(current_geo(1), current_geo(2), ...
            current_geo(3), midpoints);
        param_curr(kreis_curr) = current_log_rho(2);
        
        % Forward operation previous value i-1
        fwd_currr = app_dc.fwd.solve(sol, 1./exp(param_curr), solver_type, dc_info); % rhoa für current_x
        fwd_currr_log = log(fwd_currr);
        
        % Likelihood previous value i-1
        switch noise_type
            case 1
                % Additive Log-Like
                log_like_curr = - (1/(2*sigma^2)) * ...
                    sum((rhoa_n - fwd_currr_log).^2);
            case 2
                % Relative Log-Like
                log_like_curr = - sum(log(fwd_currr_log)) - (1/(2*sigma^2)) * ...
                    sum((rhoa_n_log - fwd_currr_log).^2 ./ fwd_currr_log.^2);
            case 3
                % Superposition
                log_like_curr = - sum(log(fwd_currr_log * noise + noise_geo')) - ...
                    0.5 * sum((rhoa_n_log - fwd_currr_log).^2 ./ ...
                    (fwd_currr_log.^2 * noise^2 + noise_geo'.^2));
        end 
    end

    % Posterior density
    post_dens = sum(log(prior_prop(:))) + log_like_prop;

    % Probability of move
    A = exp(post_dens - sum(log(prior_curr)) - log_like_curr);

    % Generate uniformly distributed random variable
    u = rand(1,1);
    
    if u <= min(A, 1)
        % Random walk
        target_leitf(:,i) = 1./exp(proposed_log_rho);
        target_geo(:,i) = proposed_geo;
        
        % Results
        result.layer(i-1,:) = exp(proposed_log_rho);
        result.geo(i-1,:)= proposed_geo;
        result.stat(i-1,1) = log_like_prop;
        result.stat(i-1,2) = post_dens;
        
        % Auxiliary variable
        mittelA = mittelA + 1;
        hlp_time = 1;
    else
        % Random walk
        target_leitf(:,i) = 1./exp(current_log_rho); 
        target_geo(:,i) = current_geo;
        
        % Auxiliary variable
        hlp_time = 0;
    end
    if i == round(0.25*n)
        fprintf('Progress: 25 %% \n');
    end
    if i == round(0.5*n)
        fprintf('Progress: 50 %% \n');
    end
    if i == round(0.75*n)
        fprintf('Progress: 75 %% \n');
    end
end

% Acceptance probability of random walk
akwrs = (mittelA/n) * 100;

%% Estimator 1 - MAP

% Sort for maximal posteriori density
[~, s_idx] = sort(result.stat(:, end), 'descend'); 
result.stat = result.stat(s_idx, :); 
result.layer = result.layer(s_idx, :);
result.geo = result.geo(s_idx, :);

idx = find(result.stat(:,1), 1); 
result.stat = result.stat(idx:end, :); 
result.layer = result.layer(idx:end, :);
result.geo = result.geo(idx:end, :);

% Computation of the map_n most probable modell answers
map_n = 2;
model_map = zeros(len_split, map_n);

for i = 1:map_n
    model_map(:,i) = result.layer(i,:);
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
end_n = n;
model_mean2 = mean(1./target_leitf(:,burn:end_n),2);
model_var2 = var(1./target_leitf(:,burn:end_n),0,2);
model_std2 = std(1./target_leitf(:,burn:end_n),0,2);
model_stat2 = [model_mean2 model_std2];

x_mean = mean(target_geo(1,burn:end_n));
y_mean = mean(target_geo(2,burn:end_n));
r_mean = mean(target_geo(3,burn:end_n));
geo_mean = [x_mean; y_mean; r_mean];
geo_std = std(target_geo(:,burn:end_n),0,2);
geo_stat = [geo_mean geo_std];

% Quantile 
rho_quant2 = quantile(1./target_leitf(:,burn:end_n), ...
    [0.025 0.25 0.50 0.75 0.975],2);

%% Random Walk Plot

rho = [800, 300];
rho_c = {'Hintergrund', 'Anomalie'};
geo = [2, -13,5];
geo_c = ['x', 'z', 'r'];

maxxil = max(1./target_leitf, [], 'all');
minnil = min(1./target_leitf, [], 'all');

figure(3)
for m = 1:len_split
    subplot(2,3,m)
        hold on
            plot(1:n, 1./target_leitf(m,:), 'ko'); 
            plot(1:n, ones(1,n)*rho(1,m), 'LineWidth', 2.5);
            plot([burn, burn], [minnil, maxxil], '-b', 'LineWidth', 2.5);
        set(gca, 'Fontsize', 15)
        set(gca, 'YScale', 'log')
        title(rho_c(1,m) +" (" + rho(1,m) + "\Omegam)")
        xlabel('Index')
        if m == 1
            ylabel('Zustand \rho [\Omegam]')
            legend('Zustand', 'Wahrer Wert', 'burn in', 'location', 'southeast')
        end
        xlim([0,n])
        xticks(0:10000:40000)
        ylim([minnil, maxxil])
        grid on
end
    subplot(2,3,4)
    m = 3;
    hold on
        plot(1:n, target_geo(m-len_split,:), 'ko'); 
        plot(1:n, ones(1,n)*geo(1,m-len_split), 'LineWidth', 2.5);
        plot([burn, burn], [-18 18], '-b', 'LineWidth', 2.5);
    hold off
    set(gca, 'Fontsize', 15)
    title("" + geo_c(m-len_split) + " (" + ...
        round(geo(1),1) + "m)")
    xticks(0:10000:40000)
    xlabel('Index')
    xlim([0,n])
    ylabel('Zustand \rho [\Omegam]')
    ylim([-18 , 18])
    yticks([-18, -10, 0, 10, 18])
    grid minor
    
    subplot(2,3,5)
    m = 4;
    hold on
        plot(1:n, target_geo(m-len_split,:), 'ko'); 
        plot(1:n, ones(1,n)*geo(1,2), 'LineWidth', 2.5);
        plot([burn, burn], [-18 18], '-b', 'LineWidth', 2.5);
    hold off
    set(gca, 'Fontsize', 15)
    title("" + geo_c(m-len_split) + " (" + ...
        round(geo(m-len_split),1) + "m)")
    xticks(0:10000:40000)
    xlabel('Index')
    xlim([0,n])
    ylabel('Zustand \rho [\Omegam]')
    ylim([-18 , -12])
    yticks([-18, -15, -12])
    grid minor
    
    subplot(2,3,6)
    m = 5;
    hold on
        plot(1:n, target_geo(m-len_split,:), 'ko'); 
        plot([burn, burn], [-18 18], '-b', 'LineWidth', 2.5);
        plot(1:n, ones(1,n)*geo(1,m-len_split), '-r', 'LineWidth', 2.5);
    hold off
    set(gca, 'Fontsize', 15)
    title("" + geo_c(m-len_split) + " (" + ...
        round(geo(m-len_split),1) + "m)")
    xticks(0:10000:40000)
    xlabel('Index')
    xlim([0,n])
    ylabel('Zustand \rho [\Omegam]')
    ylim([4, 10])
    yticks([4, 6, 8, 10])
    grid minor

%% Priori-distribution & ksdensity

prio_hist = zeros(1, n);
for i = 1:n
    prio_hist(i) = lognrnd(5.5, 0.9);
end

%% Comparison Priori & Posteriori

for k = 1:len_split
    [f(k,:), xi(k,:)] = ksdensity(1./target_leitf(k,:), (1:2:4000));
end
for k = 1:len_geo
    [fg(k,:), xig(k,:)] = ksdensity(target_geo(k,:), -18:0.3:20);
end
prior_histks = ksdensity(prio_hist, (1:2:4000));
xj = (-20:0.5:12);

xo = 18;
xu = -18;
x = rand(1)*(xo-xu)+xu;
zo = -12;
zu = -18;
z = rand(1)*(zo-zu)+zu; 
r1 = 4; 
r2 = 10; 
r = abs(rand(1)*(r1-r2)+r2);

for i = 1:2*n
    prio_hist(i) = lognrnd(5.5, 0.9); 
    prio_histx(i) = rand(1)*(xo-xu)+xu; 
    prio_histz(i) = rand(1)*(zo-zu)+zu;
    prio_histr(i) = rand(1)*(4-10)+10;
end


mp = [800, 300, 4, 10, 8];

xix = (-18:1:18);
xiz = (-18:0.3:-12);
xir = (4:0.3:10);

[fgx, xgx] = ksdensity(prio_histx, xix);
[fgz, xgz] = ksdensity(prio_histz, xiz);
[fgr, xgr] = ksdensity(prio_histr, xir);

%% Comparison Priori & Posteriori Plot

figure(4)
subplot(2,3,1)
set(gca, 'Fontsize', 15)
hold on
plot(xi(1,:), f(1,:), 'LineWidth', 2.5)
plot(xi(1,:), prior_histks, 'LineWidth', 2.5)
set(gca, 'XScale', 'log')
title('Half space 800 \Omegam')
xlabel('Spez. el. Widerstand [ \Omegam]')
legend("Posteriori", "Priori", 'location', 'northwest', 'Fontsize', 11)
xlim([10 4000])
xticks(logspace(1,4,4))
grid on
ax = gca;
ax.YAxis.Exponent = 0;
    
subplot(2,3,2)
set(gca, 'Fontsize', 15)
hold on
plot(xi(1,:), f(2,:), 'LineWidth', 2.5)
plot(xi(1,:), prior_histks, 'LineWidth', 2.5)
set(gca, 'XScale', 'log')
title('Anomaly 300 \Omegam')
xlabel('Resistivity [\Omegam]')
xlim([10 4000])
xticks(logspace(1,4,4))
grid on
ax = gca;
ax.YAxis.Exponent = 0;

subplot(2,3,4)
set(gca, 'Fontsize', 15)
hold on
plot(xig(1,:), fg(1,:), 'LineWidth', 2.5)
plot(xgx, fgx, 'LineWidth', 2.5)
xlim([-18, 18])
title("x " + round(geo(1),1) + " m")
xlabel('Profile position [m]')
grid on

subplot(2,3,5)
set(gca, 'Fontsize', 15)
hold on
plot(xig(2,:), fg(2,:), 'LineWidth', 2.5)
plot(xgz, fgz, 'LineWidth', 2.5)
xlim([-18, -12])
title("z " + round(geo(2),1) + " m")
xlabel('Depth [m]')
grid on

subplot(2,3,6)
set(gca, 'Fontsize', 15)
hold on
plot(xig(3,:), fg(3,:), 'LineWidth', 2.5)
plot(xgr, fgr, 'LineWidth', 2.5)
xlim([4, 10])
xticks([4,6,8,10])
title("r " + round(geo(3),1) + " m")
xlabel('Radius [m]')
grid on

%% Error plots

param_result = ones(size(param));
param_result(:) = model_stat2(1,1);
kreis_result = anomaly_triangle(geo_stat(1,1), geo_stat(2,1), ...
    geo_stat(3,1), midpoints);
param_result(kreis_result) = model_stat2(2,1);

hlp_result = app_dc.fwd.solve(sol, 1./param_result, solver_type, dc_info);
err_rhoa = abs(hlp_result - rhoa);

% MAP
param_result_map = ones(size(param));
param_result_map(:) = result.layer(1,1);
kreis_result_map = anomaly_triangle(result.geo(1,1), result.geo(1,2), ...
    result.geo(1,3), midpoints);
param_result_map(kreis_result_map) = result.layer(1,2);

hlp_result_map = app_dc.fwd.solve(sol, 1./param_result_map, solver_type, dc_info);
err_rhoa_map = abs(hlp_result_map - rhoa);

figure(5)
app_dc.data.plot_pseudosection(dc_info.survey, rhoa)
title('Messwert \rho_a')
set(gca,'ColorScale','log')
set(gca, 'Fontsize', 16)
setcolorbartitle('\rho_a [\Omegam]')

figure(6)
app_dc.data.plot_pseudosection(dc_info.survey, hlp_result)
title('Results mean')
% set(gca,'ColorScale','log')
set(gca, 'Fontsize', 16)
setcolorbartitle('\rho_a [\Omegam]')

figure(7)
app_dc.data.plot_pseudosection(dc_info.survey, err_rhoa)
title('Absoluter Fehler')
set(gca, 'Fontsize', 16)
setcolorbartitle('\rho_a [\Omegam]')

%% Plot model outcomes

paramr = 1./param;
comin = min([paramr; param_result], [], 'all');
comax = max([paramr; param_result], [], 'all');

figure(8)
plot_model(mesh, param, 'Model', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
caxis([comin comax])
setcolorbar('\rho [\Omegam]', comin, comax)
xlim([-25 25])
ylim([-25 2])

figure(9)
plot_model(mesh, 1./param_result, 'Result mean', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
setcolorbar('\rho [\Omegam]', comin, comax)
caxis([comin comax])
xlim([-25 25])
ylim([-25 2])

figure(10)
plot_model(mesh, 1./param_result_map, 'Resultat MAP', dc_info)
set(gca, 'Fontsize', 16)
set(gca,'ColorScale','log')
setcolorbar('\rho [\Omegam]', comin, comax)
caxis([comin comax])
xlim([-25 25])
ylim([-25 2])

%% Anomaly with often used triangles

for t = 1:(n-burn-1)
    kreisbol(t,:) = anomaly_triangle(target_geo(1,burn+t), ...
        target_geo(2,burn+t), ...
        target_geo(3,burn+t), midpoints);
end

kreisbol_sum = sum(kreisbol,1);
kreisbol_grey = kreisbol_sum / max(kreisbol_sum);
kreisbol_grey = kreisbol_grey';

figure(11)
plot_model_gray(mesh, 1./kreisbol_grey, 'Oft verwendete Elemente', dc_info)
set(gca, 'Fontsize', 16)
hcb = colorbar;
colorTitleHandle = get(hcb,'Title');
set(colorTitleHandle ,'String','Prozent');
xlim([-25 25])
ylim([-25 2])

%% Fixed depth between 13-14 and 15-16
hlp_val = target_geo(2,burn:end);
B = hlp_val(hlp_val < -13 & hlp_val > -14);
C = find(hlp_val < -13 & hlp_val > -14);
hlp_val2 = 1./target_leitf(2,burn:end);
hlp_val3 = target_geo(3,burn:end);
rhoano = hlp_val2(C);
flaeche = hlp_val3(C);

B2 = hlp_val(hlp_val < -15 & hlp_val > -16);
C2 = find(hlp_val < -15 & hlp_val > -16);
hlp_val4 = 1./target_leitf(2,burn:end);
hlp_val5 = target_geo(3,burn:end);
rhoano2 = hlp_val4(C2);
flaeche2 = hlp_val5(C2);

figure(12)
scatter(flaeche.^2, rhoano, 'b.')
hold on
set(gca, 'Fontsize', 16)
title('Crossplot r^2 zu \rho_{ano}, festgehaltene Tiefe zw. 15 und 16')
xlabel('Fläche r^2 [m^2]')
ylabel('\rho [\Omegam]')
grid on

%% Fixed radius
hlp_val2 = target_geo(3,burn:end);
hlp_val6 = 1./target_leitf(2,burn:end);
hlp_val7= target_geo(2,burn:end);

C45 = find(hlp_val2 > 4 & hlp_val2 < 5);
rhoano45 = hlp_val6(C45);
tiefe45 = abs(hlp_val7(C45));

C56 = find(hlp_val2 > 5 & hlp_val2 < 6);
rhoano56 = hlp_val6(C56);
tiefe56 = abs(hlp_val7(C56));

C67 = find(hlp_val2 > 6 & hlp_val2 < 7);
rhoano67 = hlp_val6(C67);
tiefe67 = abs(hlp_val7(C67));

C78 = find(hlp_val2 > 7 & hlp_val2 < 8);
rhoano78 = hlp_val6(C78);
tiefe78 = abs(hlp_val7(C78));

C89 = find(hlp_val2 > 8 & hlp_val2 < 9);
rhoano89 = hlp_val6(C89);
tiefe89 = abs(hlp_val7(C89));

C91 = find(hlp_val2 > 9 & hlp_val2 < 10);
rhoano91 = hlp_val6(C91);
tiefe91 = abs(hlp_val7(C91));

figure(64)
scatter(tiefe45, rhoano45, '.')
hold on
scatter(tiefe56, rhoano56, '.')
scatter(tiefe67, rhoano67, '.')
scatter(tiefe78, rhoano78, '.')
scatter(tiefe89, rhoano89, '.')
scatter(tiefe91, rhoano91, '.')
set(gca, 'Fontsize', 16)
title('Streudiagramm Tiefe zu \rho_{ano} mit festem Radius')
xlabel('Tiefe z [m]')
ylabel('\rho_{ano} [\Omegam]')
legend('4-5 m','5-6 m','6-7 m','7-8 m','8-9 m','9-10 m','location','best')
grid on

%% Crossplots
K1 = [ones(length(target_geo(3,burn:end)),1) (pi*target_geo(3,burn:end).^2)'];
K2 = ( 1./target_leitf(2,burn:end))';
K = K1 \K2;
yCalc1 = K1*K;

figure(13)
scatter(pi*target_geo(3,burn:end).^2, 1./target_leitf(2,burn:end), '.')
hold on
scatter(pi*result.geo(1:1000,3).^2, result.layer(1:1000,2), 'r.')
plot(pi*target_geo(3,burn:end).^2, yCalc1, 'k-', 'Linewidth', 2.5)
set(gca, 'Fontsize', 16)
title('Streudiagramm r^2 zu \rho_{ano}')
xlabel('Fläche r^2 [m^2]')
ylabel('\rho_{ano} [\Omegam]')
legend('Gesamter Random Walk','1000 besten MAP','Regressionsgerade', 'location', 'northwest')
grid on

figure(14)
scatter(abs(target_geo(2,burn:end)), 1./target_leitf(2,burn:end), '.')
set(gca, 'Fontsize', 16)
title('Crossplot z zu \rho_{ano}')
xlabel('Tiefe z [m]')
ylabel('\rho [\Omegam]')
grid on

figure(15)
scatter(abs(target_geo(2,burn:end)), pi*target_geo(3,burn:end).^2, '.')
set(gca, 'Fontsize', 16)
title('Crossplot Tiefe zu Fläche')
xlabel('Tiefe z [m]')
ylabel('Fläche r^2 [m^2]')
grid on

%% Time display
et = toc;
etchar = sprintf('Elapsed time: %.2f sec. \n', et);
fprintf(etchar);

%% Save option

% save('variables_zylinderanomaly.mat')

%% Helpers

function plot_model(mesh, sigma, str, dc_info)

    meshing.plot_mesh(mesh, 'cell_markers', 1./sigma);
    ylabel('z in m');
    title(str);
    xlabel('x in m');
    hold on
        plot(dc_info.survey.ele(:, 1), dc_info.survey.ele(:, 2), 'ro');
        meshing.plot_mesh(mesh)
    hold off
    drawnow;
end

function kreis = anomaly_triangle(x, y, r, midpoints)
    cycle = [x; y]; 
    middl = vecnorm(cycle - midpoints);
    kreis = middl <= r; 
end

function plot_model_gray(mesh, sigma, str, dc_info)

    meshing.plot_mesh(mesh, 'cell_markers', 1./sigma);
    ylabel('z [m]');
    title(str);
    xlabel('x [m]');
    hold on
        plot(dc_info.survey.ele(:, 1), dc_info.survey.ele(:, 2), 'ro');
        meshing.plot_mesh(mesh)
    hold off
    drawnow;
    colormap(flipud(gray))
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
    set(hcb,'YTick', round(logspace(log10(comin), log10(comax),5),-2))
end