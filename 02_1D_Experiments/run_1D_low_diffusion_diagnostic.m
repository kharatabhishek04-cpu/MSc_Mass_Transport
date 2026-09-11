% 1D low-diffusion stability diagnostic
clear; clc; close all;

% Load benchmark functions
project_folder = fileparts(fileparts(mfilename('fullpath')));
benchmark_folder = fullfile(project_folder,'01_1D_Benchmark');
addpath(benchmark_folder);

fprintf('LBM functions loaded from: %s\n\n', benchmark_folder);

% Domain
channel_length = 400;
Lx = 401;
Ly = 5;

dx = channel_length/(Lx-1);
dy = dx;
dt = 0.1;
e = dx/dt;

% Model parameters
tau = 1.1;
u0 = 0.01;
v0 = 0.0;

C0 = 1.0;
x0 = 10.0;
source_index = round(x0/dx) + 1;

% Diffusion cases
diffusion_values = [0.006; 0.005; 0.004; 0.003; 0.002];
number_of_cases = length(diffusion_values);

% Times used for checking the solution
snapshot_times = [100; 1000; 5000; 10000];
snapshot_iterations = round(snapshot_times/dt);
number_of_snapshots = length(snapshot_times);

final_time = max(snapshot_times);
max_iterations = round(final_time/dt);

negative_tolerance = -1e-12;
blowup_limit = 10;

% Estimated positivity threshold
D_critical = u0*(tau - 0.5)*dx;

fprintf('LOW-DIFFUSION STABILITY DIAGNOSTIC\n');
fprintf('Velocity = %.4f m/s\n', u0);
fprintf('Tau = %.2f\n', tau);
fprintf('dx = %.4f m\n', dx);
fprintf('dt = %.4f s\n', dt);
fprintf('Estimated D threshold = %.6f m^2/s\n\n', D_critical);

fprintf('Diffusion values:\n');
disp(diffusion_values');

fprintf('Diagnostic times:\n');
disp(snapshot_times');

% Coordinates
x_coordinates = (0:Lx-1)'*dx;
middle_y = ceil(Ly/2);

% Case summary storage
Diffusion = zeros(number_of_cases,1);
Peclet = zeros(number_of_cases,1);
InitialMinimumFeq = zeros(number_of_cases,1);
InitialNegativePopulationCount = zeros(number_of_cases,1);
FirstNegativeTime_s = NaN(number_of_cases,1);
MostNegativeConcentration = zeros(number_of_cases,1);
LargestAbsoluteConcentration = zeros(number_of_cases,1);
FinalReachedTime_s = zeros(number_of_cases,1);
FinalMass = zeros(number_of_cases,1);
StoppedEarly = false(number_of_cases,1);
Status = strings(number_of_cases,1);

% Snapshot results
maximum_snapshot_rows = number_of_cases*number_of_snapshots;

Snapshot_Diffusion = zeros(maximum_snapshot_rows,1);
Snapshot_Peclet = zeros(maximum_snapshot_rows,1);
Snapshot_Time_s = zeros(maximum_snapshot_rows,1);
Snapshot_MinimumC = zeros(maximum_snapshot_rows,1);
Snapshot_MaximumC = zeros(maximum_snapshot_rows,1);
Snapshot_RMSE = zeros(maximum_snapshot_rows,1);
Snapshot_MAE = zeros(maximum_snapshot_rows,1);
Snapshot_R2 = zeros(maximum_snapshot_rows,1);
Snapshot_Mass = zeros(maximum_snapshot_rows,1);
Snapshot_NegativeCount = zeros(maximum_snapshot_rows,1);
Snapshot_HasNegativeC = false(maximum_snapshot_rows,1);

snapshot_row = 0;

% Save min/max concentration every 10 seconds
history_interval = 100;
history_points = floor(max_iterations/history_interval);

history_time = NaN(history_points,number_of_cases);
history_minC = NaN(history_points,number_of_cases);
history_maxC = NaN(history_points,number_of_cases);

% Final profiles
final_profiles = NaN(Lx,number_of_cases);
final_profile_times = zeros(number_of_cases,1);

for d_index = 1:number_of_cases

    D = diffusion_values(d_index);
    Pe = u0*dx/D;

    fprintf('\nTesting D = %.4f m^2/s\n', D);
    fprintf('Pe = %.4f\n', Pe);

    comparison_tolerance = 1e-12;

    if D < D_critical - comparison_tolerance
        fprintf('Warning: D is below estimated positivity threshold.\n');
    else
        fprintf('D is at or above estimated positivity threshold.\n');
    end

    % Reset model for each diffusion value
    solid = zeros(Lx,Ly);
    Cen = zeros(Lx,Ly);

    u = u0*ones(Lx,Ly);
    v = v0*ones(Lx,Ly);

    Cen(source_index,:) = C0;

    [ex,ey] = setup(e);

    feq = compute_feq( ...
        u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

    initial_min_feq = min(feq(:));
    initial_negative_count = sum(feq(:) < negative_tolerance);

    fprintf('Initial minimum feq = %.12e\n', initial_min_feq);
    fprintf('Initial negative populations = %d\n', initial_negative_count);

    f = feq;

    first_negative_time = NaN;
    most_negative_C = min(Cen(:));
    largest_absolute_C = max(abs(Cen(:)));

    stopped_early = false;
    status_text = "Completed";

    snapshot_number = 1;
    history_index = 0;
    last_iteration = 0;

    for iteration = 1:max_iterations

        ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);

        ftemp = Peri_BC_inOut(2,Ly,ftemp,Lx);
        ftemp = BC_NorthSouth(ftemp,1,Lx,1,Ly);

        % Copy lower boundary
        for x = 1:Lx
            for a = 1:5
                ftemp(a,x,1) = ftemp(a,x,2);
            end
        end

        % Copy upper boundary
        for x = 1:Lx
            for a = 1:5
                ftemp(a,x,Ly) = ftemp(a,x,Ly-1);
            end
        end

        [Cen,f] = solution(ftemp,1,Lx,1,Ly);

        feq = compute_feq( ...
            u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

        current_time = iteration*dt;
        last_iteration = iteration;

        current_min_C = min(Cen(:));
        current_max_abs_C = max(abs(Cen(:)));

        most_negative_C = min(most_negative_C,current_min_C);
        largest_absolute_C = max(largest_absolute_C,current_max_abs_C);

        % Record the first negative concentration
        if isnan(first_negative_time) && current_min_C < negative_tolerance

            first_negative_time = current_time;

            fprintf('\nFirst negative concentration detected\n');
            fprintf('D = %.4f\n', D);
            fprintf('Time = %.2f s\n', current_time);
            fprintf('Min C = %.12e\n\n', current_min_C);
        end

        % Store concentration history
        if mod(iteration,history_interval) == 0

            history_index = history_index + 1;

            history_time(history_index,d_index) = current_time;
            history_minC(history_index,d_index) = current_min_C;
            history_maxC(history_index,d_index) = max(Cen(:));
        end

        % Snapshot calculations
        if snapshot_number <= number_of_snapshots && ...
                iteration == snapshot_iterations(snapshot_number)

            numerical_profile = Cen(:,middle_y);

            analytical_profile = ...
                C0./sqrt(4*pi*D*current_time) .* ...
                exp(-((x_coordinates-x0-u0*current_time).^2) ./ ...
                (4*D*current_time));

            error_values = numerical_profile - analytical_profile;

            current_RMSE = sqrt(mean(error_values.^2));
            current_MAE = mean(abs(error_values));

            SS_res = sum(error_values.^2);
            SS_tot = sum((analytical_profile - mean(analytical_profile)).^2);

            current_R2 = 1 - SS_res/SS_tot;

            current_mass = sum(numerical_profile)*dx;

            current_negative_count = ...
                sum(numerical_profile < negative_tolerance);

            current_negative_flag = current_negative_count > 0;

            snapshot_row = snapshot_row + 1;

            Snapshot_Diffusion(snapshot_row) = D;
            Snapshot_Peclet(snapshot_row) = Pe;
            Snapshot_Time_s(snapshot_row) = current_time;
            Snapshot_MinimumC(snapshot_row) = min(numerical_profile);
            Snapshot_MaximumC(snapshot_row) = max(numerical_profile);
            Snapshot_RMSE(snapshot_row) = current_RMSE;
            Snapshot_MAE(snapshot_row) = current_MAE;
            Snapshot_R2(snapshot_row) = current_R2;
            Snapshot_Mass(snapshot_row) = current_mass;
            Snapshot_NegativeCount(snapshot_row) = current_negative_count;
            Snapshot_HasNegativeC(snapshot_row) = current_negative_flag;

            fprintf('\nD = %.4f | Time = %.0f s\n', D,current_time);
            fprintf('Min C = %.12e\n', min(numerical_profile));
            fprintf('Max C = %.10f\n', max(numerical_profile));
            fprintf('Negative points = %d\n', current_negative_count);
            fprintf('RMSE = %.10f\n', current_RMSE);
            fprintf('MAE = %.10f\n', current_MAE);
            fprintf('R2 = %.10f\n', current_R2);
            fprintf('Mass = %.10f\n', current_mass);

            snapshot_number = snapshot_number + 1;
        end

        % Stop if simulation fails
        if any(~isfinite(Cen(:)))

            fprintf('\nSimulation stopped: NaN or Inf detected.\n');

            stopped_early = true;
            status_text = "Stopped: NaN/Inf";
            break;
        end

        if current_max_abs_C > blowup_limit

            fprintf('\nSimulation stopped: |C| exceeded %.2f.\n', ...
                blowup_limit);

            stopped_early = true;
            status_text = "Stopped: severe numerical blow-up";
            break;
        end
    end

    % Final result for this D
    final_profiles(:,d_index) = Cen(:,middle_y);

    final_reached_time = last_iteration*dt;
    final_profile_times(d_index) = final_reached_time;

    final_mass = sum(Cen(:,middle_y))*dx;

    Diffusion(d_index) = D;
    Peclet(d_index) = Pe;
    InitialMinimumFeq(d_index) = initial_min_feq;
    InitialNegativePopulationCount(d_index) = initial_negative_count;
    FirstNegativeTime_s(d_index) = first_negative_time;
    MostNegativeConcentration(d_index) = most_negative_C;
    LargestAbsoluteConcentration(d_index) = largest_absolute_C;
    FinalReachedTime_s(d_index) = final_reached_time;
    FinalMass(d_index) = final_mass;
    StoppedEarly(d_index) = stopped_early;
    Status(d_index) = status_text;

    fprintf('\nCase complete: D = %.4f\n', D);
    fprintf('Initial minimum feq = %.12e\n', initial_min_feq);
    fprintf('Initial negative f count = %d\n', initial_negative_count);

    if isnan(first_negative_time)
        fprintf('First negative C time = NONE\n');
    else
        fprintf('First negative C time = %.2f s\n', first_negative_time);
    end

    fprintf('Most negative C = %.12e\n', most_negative_C);
    fprintf('Largest |C| = %.10f\n', largest_absolute_C);
    fprintf('Final reached time = %.2f s\n', final_reached_time);
    fprintf('Final mass = %.10f\n', final_mass);
    fprintf('Status = %s\n', status_text);
end

% Remove unused snapshot rows
Snapshot_Diffusion = Snapshot_Diffusion(1:snapshot_row);
Snapshot_Peclet = Snapshot_Peclet(1:snapshot_row);
Snapshot_Time_s = Snapshot_Time_s(1:snapshot_row);
Snapshot_MinimumC = Snapshot_MinimumC(1:snapshot_row);
Snapshot_MaximumC = Snapshot_MaximumC(1:snapshot_row);
Snapshot_RMSE = Snapshot_RMSE(1:snapshot_row);
Snapshot_MAE = Snapshot_MAE(1:snapshot_row);
Snapshot_R2 = Snapshot_R2(1:snapshot_row);
Snapshot_Mass = Snapshot_Mass(1:snapshot_row);
Snapshot_NegativeCount = Snapshot_NegativeCount(1:snapshot_row);
Snapshot_HasNegativeC = Snapshot_HasNegativeC(1:snapshot_row);

% Summary tables
case_summary = table( ...
    Diffusion, ...
    Peclet, ...
    InitialMinimumFeq, ...
    InitialNegativePopulationCount, ...
    FirstNegativeTime_s, ...
    MostNegativeConcentration, ...
    LargestAbsoluteConcentration, ...
    FinalReachedTime_s, ...
    FinalMass, ...
    StoppedEarly, ...
    Status);

snapshot_summary = table( ...
    Snapshot_Diffusion, ...
    Snapshot_Peclet, ...
    Snapshot_Time_s, ...
    Snapshot_MinimumC, ...
    Snapshot_MaximumC, ...
    Snapshot_NegativeCount, ...
    Snapshot_RMSE, ...
    Snapshot_MAE, ...
    Snapshot_R2, ...
    Snapshot_Mass, ...
    Snapshot_HasNegativeC);

fprintf('\nLOW-DIFFUSION CASE SUMMARY\n');
disp(case_summary);

fprintf('\nSNAPSHOT SUMMARY\n');
disp(snapshot_summary);

% Save results
results_directory = fullfile( ...
    project_folder,'06_Results','1D_Low_Diffusion_Diagnostic');

if ~exist(results_directory,'dir')
    mkdir(results_directory);
end

writetable(case_summary, ...
    fullfile(results_directory,'1D_low_diffusion_case_summary.csv'));

writetable(snapshot_summary, ...
    fullfile(results_directory,'1D_low_diffusion_snapshot_summary.csv'));

save(fullfile(results_directory,'1D_low_diffusion_diagnostic.mat'));

% Minimum concentration through time
figure;
hold on;

for d_index = 1:number_of_cases

    valid_history = ~isnan(history_time(:,d_index));

    plot( ...
        history_time(valid_history,d_index), ...
        history_minC(valid_history,d_index), ...
        'LineWidth',1.5);
end

yline(negative_tolerance,'--','Negative threshold');

xlabel('Time (s)');
ylabel('Minimum concentration');
title('Low-Diffusion Stability Diagnostic');

legend( ...
    'D=0.006', ...
    'D=0.005', ...
    'D=0.004', ...
    'D=0.003', ...
    'D=0.002', ...
    'Location','best');

grid on;
hold off;

exportgraphics( ...
    gcf, ...
    fullfile(results_directory, ...
    '1D_low_diffusion_minimum_concentration.png'), ...
    'Resolution',300);

% Final concentration profiles
figure;
hold on;

for d_index = 1:number_of_cases
    plot(x_coordinates,final_profiles(:,d_index),'LineWidth',1.5);
end

xlabel('x (m)');
ylabel('C (kg/m)');
title('Low-Diffusion Final Concentration Profiles');

legend( ...
    'D=0.006', ...
    'D=0.005', ...
    'D=0.004', ...
    'D=0.003', ...
    'D=0.002', ...
    'Location','best');

grid on;
hold off;

exportgraphics( ...
    gcf, ...
    fullfile(results_directory, ...
    '1D_low_diffusion_final_profiles.png'), ...
    'Resolution',300);

fprintf('\nLow-diffusion diagnostic completed.\n');
fprintf('Estimated D threshold = %.6f m^2/s\n',D_critical);
fprintf('Results saved in: %s\n',results_directory);