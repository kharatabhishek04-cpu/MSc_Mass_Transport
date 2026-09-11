% 1D velocity experiment
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

tau = 1.1;

% Fixed diffusion
D = 0.01;
v0 = 0.0;

% Source
C0 = 1.0;
x0 = 10.0;
source_index = round(x0/dx) + 1;

% Velocity cases
velocity_values = [0.0050; 0.0075; 0.0100; 0.0125; 0.0150];
number_of_velocities = length(velocity_values);

% Times to record results
snapshot_times = [5000; 10000; 15000; 20000];
snapshot_iterations = round(snapshot_times/dt);
number_of_snapshots = length(snapshot_times);

final_time = max(snapshot_times);
max_iterations = round(final_time/dt);

% Approximate positivity limit
u_critical = D/((tau - 0.5)*dx);

fprintf('1D CONTROLLED VELOCITY EXPERIMENT\n');
fprintf('D = %.4f m^2/s\n', D);
fprintf('tau = %.2f\n', tau);
fprintf('dx = %.4f m\n', dx);
fprintf('dt = %.4f s\n', dt);
fprintf('Estimated velocity limit = %.6f m/s\n\n', u_critical);

fprintf('Velocity values: ');
fprintf('%.4f ', velocity_values);
fprintf('\n');

fprintf('Observation times: ');
fprintf('%.0f ', snapshot_times);
fprintf('s\n\n');

xc = (0:Lx-1)'*dx;
mid_y = ceil(Ly/2);

% Result storage
total_cases = number_of_velocities*number_of_snapshots;

CaseID = strings(total_cases,1);
Velocity = zeros(total_cases,1);
Peclet = zeros(total_cases,1);
Time_s = zeros(total_cases,1);

RMSE = zeros(total_cases,1);
MAE = zeros(total_cases,1);
R2 = zeros(total_cases,1);
MaxAbsoluteError = zeros(total_cases,1);
MinimumConcentration = zeros(total_cases,1);

NumericalPeak = zeros(total_cases,1);
AnalyticalPeak = zeros(total_cases,1);
NumericalPeakX_m = zeros(total_cases,1);
AnalyticalPeakX_m = zeros(total_cases,1);

NumericalCentroidX_m = zeros(total_cases,1);
ExpectedCentroidX_m = zeros(total_cases,1);

NumericalSigma_m = zeros(total_cases,1);
TheoreticalSigma_m = zeros(total_cases,1);

NumericalMass = zeros(total_cases,1);
AnalyticalMass = zeros(total_cases,1);

HasNegativeConcentration = false(total_cases,1);
Runtime_seconds = zeros(total_cases,1);

% Profiles for CSV export
profile_velocity = [];
profile_Pe = [];
profile_time = [];
profile_x = [];
profile_numerical = [];
profile_analytical = [];
profile_error = [];

% Profiles for plots
all_numerical_profiles = zeros( ...
    Lx,number_of_velocities,number_of_snapshots);

all_analytical_profiles = zeros( ...
    Lx,number_of_velocities,number_of_snapshots);

row = 0;

for u_index = 1:number_of_velocities

    u0 = velocity_values(u_index);
    Pe = u0*dx/D;

    fprintf('\nVelocity case u = %.4f m/s (Pe = %.4f)\n', u0,Pe);

    if u0 > u_critical
        fprintf('Warning: velocity exceeds estimated positivity limit.\n');
    else
        fprintf('Velocity is below estimated positivity limit.\n');
    end

    % Reset model
    solid = zeros(Lx,Ly);
    Cen = zeros(Lx,Ly);

    u = u0*ones(Lx,Ly);
    v = v0*ones(Lx,Ly);

    Cen(source_index,:) = C0;

    [ex,ey] = setup(e);

    feq = compute_feq( ...
        u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

    initial_min_feq = min(feq(:));
    initial_negative_feq = sum(feq(:) < -1e-12);

    fprintf('Initial minimum feq = %.12e\n', initial_min_feq);
    fprintf('Initial negative populations = %d\n', initial_negative_feq);

    f = feq;

    snapshot_number = 1;
    case_timer = tic;

    for iteration = 1:max_iterations

        ftemp = collide_stream(1,Lx,1,Ly,f,feq,tau);

        ftemp = Peri_BC_inOut(2,Ly,ftemp,Lx);
        ftemp = BC_NorthSouth(ftemp,1,Lx,1,Ly);

        % Copy lower and upper y boundaries
        for x = 1:Lx
            for a = 1:5
                ftemp(a,x,1) = ftemp(a,x,2);
            end
        end

        for x = 1:Lx
            for a = 1:5
                ftemp(a,x,Ly) = ftemp(a,x,Ly-1);
            end
        end

        [Cen,f] = solution(ftemp,1,Lx,1,Ly);

        feq = compute_feq( ...
            u,v,dt,tau,e,ex,ey,D,solid,Lx,Ly,Cen);

        % Record required times
        if snapshot_number <= number_of_snapshots && ...
                iteration == snapshot_iterations(snapshot_number)

            current_time = snapshot_times(snapshot_number);
            row = row + 1;

            Cnum = Cen(:,mid_y);

            Cana = C0./sqrt(4*pi*D*current_time) .* ...
                exp(-((xc-x0-u0*current_time).^2) ./ ...
                (4*D*current_time));

            err = Cnum - Cana;

            cur_RMSE = sqrt(mean(err.^2));
            cur_MAE = mean(abs(err));

            SS_res = sum(err.^2);
            SS_tot = sum((Cana-mean(Cana)).^2);

            cur_R2 = 1 - SS_res/SS_tot;
            cur_max_error = max(abs(err));

            cur_min_C = min(Cnum);
            negative_flag = cur_min_C < -1e-12;

            % Peaks
            [num_peak,num_peak_idx] = max(Cnum);
            [ana_peak,ana_peak_idx] = max(Cana);

            num_peak_x = xc(num_peak_idx);
            ana_peak_x = xc(ana_peak_idx);

            % Mass
            num_mass = sum(Cnum)*dx;
            ana_mass = sum(Cana)*dx;

            % Centroid
            num_centroid = ...
                sum(xc.*Cnum)*dx/num_mass;

            expected_centroid = x0 + u0*current_time;

            % Plume spread
            num_sigma = sqrt( ...
                sum(((xc-num_centroid).^2).*Cnum)*dx/num_mass);

            theoretical_sigma = sqrt(2*D*current_time);

            current_runtime = toc(case_timer);

            % Store summary
            CaseID(row) = sprintf('U%.4f_T%d',u0,current_time);

            Velocity(row) = u0;
            Peclet(row) = Pe;
            Time_s(row) = current_time;

            RMSE(row) = cur_RMSE;
            MAE(row) = cur_MAE;
            R2(row) = cur_R2;
            MaxAbsoluteError(row) = cur_max_error;
            MinimumConcentration(row) = cur_min_C;

            NumericalPeak(row) = num_peak;
            AnalyticalPeak(row) = ana_peak;

            NumericalPeakX_m(row) = num_peak_x;
            AnalyticalPeakX_m(row) = ana_peak_x;

            NumericalCentroidX_m(row) = num_centroid;
            ExpectedCentroidX_m(row) = expected_centroid;

            NumericalSigma_m(row) = num_sigma;
            TheoreticalSigma_m(row) = theoretical_sigma;

            NumericalMass(row) = num_mass;
            AnalyticalMass(row) = ana_mass;

            HasNegativeConcentration(row) = negative_flag;
            Runtime_seconds(row) = current_runtime;

            % Save profiles for plots
            all_numerical_profiles(:,u_index,snapshot_number) = Cnum;
            all_analytical_profiles(:,u_index,snapshot_number) = Cana;

            % Save profiles for CSV
            profile_velocity = [profile_velocity; repmat(u0,Lx,1)];
            profile_Pe = [profile_Pe; repmat(Pe,Lx,1)];
            profile_time = [profile_time; repmat(current_time,Lx,1)];

            profile_x = [profile_x; xc];
            profile_numerical = [profile_numerical; Cnum];
            profile_analytical = [profile_analytical; Cana];
            profile_error = [profile_error; err];

            fprintf('\nt = %.0f s\n', current_time);
            fprintf('RMSE = %.10f\n', cur_RMSE);
            fprintf('MAE = %.10f\n', cur_MAE);
            fprintf('R2 = %.10f\n', cur_R2);
            fprintf('Min C = %.12e\n', cur_min_C);
            fprintf('Numerical peak = %.10f\n', num_peak);
            fprintf('Analytical peak = %.10f\n', ana_peak);
            fprintf('Numerical x peak = %.2f m\n', num_peak_x);
            fprintf('Expected centre = %.2f m\n', expected_centroid);
            fprintf('Numerical centre = %.4f m\n', num_centroid);
            fprintf('Numerical sigma = %.4f m\n', num_sigma);
            fprintf('Theory sigma = %.4f m\n', theoretical_sigma);
            fprintf('Numerical mass = %.10f\n', num_mass);
            fprintf('Negative C = %d\n', negative_flag);

            snapshot_number = snapshot_number + 1;
        end
    end
end

% Summary table
summary_table = table( ...
    CaseID, ...
    Velocity, ...
    Peclet, ...
    Time_s, ...
    RMSE, ...
    MAE, ...
    R2, ...
    MaxAbsoluteError, ...
    MinimumConcentration, ...
    NumericalPeak, ...
    AnalyticalPeak, ...
    NumericalPeakX_m, ...
    AnalyticalPeakX_m, ...
    NumericalCentroidX_m, ...
    ExpectedCentroidX_m, ...
    NumericalSigma_m, ...
    TheoreticalSigma_m, ...
    NumericalMass, ...
    AnalyticalMass, ...
    HasNegativeConcentration, ...
    Runtime_seconds);

fprintf('\n');
disp(summary_table);

% Save results
results_directory = fullfile( ...
    project_folder,'06_Results','1D_Velocity');

if ~exist(results_directory,'dir')
    mkdir(results_directory);
end

writetable( ...
    summary_table, ...
    fullfile(results_directory,'1D_velocity_summary.csv'));

profiles_table = table( ...
    profile_velocity, ...
    profile_Pe, ...
    profile_time, ...
    profile_x, ...
    profile_numerical, ...
    profile_analytical, ...
    profile_error, ...
    'VariableNames',{ ...
    'Velocity', ...
    'Peclet', ...
    'Time_s', ...
    'x_m', ...
    'NumericalConcentration', ...
    'AnalyticalConcentration', ...
    'Error'});

writetable( ...
    profiles_table, ...
    fullfile(results_directory,'1D_velocity_profiles.csv'));

save(fullfile(results_directory,'1D_velocity_experiment.mat'));

% Effect of velocity at each time
figure;
tiledlayout(2,2);

for t_index = 1:number_of_snapshots

    nexttile;
    hold on;

    for u_index = 1:number_of_velocities
        plot( ...
            xc, ...
            all_numerical_profiles(:,u_index,t_index), ...
            'LineWidth',1.5);
    end

    xlabel('x (m)');
    ylabel('C (kg/m)');
    title(sprintf('t = %.0f s',snapshot_times(t_index)));

    xlim([0 channel_length]);
    grid on;

    if t_index == 1
        legend( ...
            'u=0.0050', ...
            'u=0.0075', ...
            'u=0.0100', ...
            'u=0.0125', ...
            'u=0.0150', ...
            'Location','best');
    end

    hold off;
end

sgtitle('1D Controlled Velocity Experiment');

exportgraphics( ...
    gcf, ...
    fullfile(results_directory,'1D_velocity_comparison.png'), ...
    'Resolution',300);

% Numerical vs analytical at 20000 s
validation_time_index = 4;

figure;
tiledlayout(2,3);

for u_index = 1:number_of_velocities

    nexttile;

    plot( ...
        xc, ...
        all_numerical_profiles(:,u_index,validation_time_index), ...
        'o','MarkerSize',3);

    hold on;

    plot( ...
        xc, ...
        all_analytical_profiles(:,u_index,validation_time_index), ...
        '-','LineWidth',1.5);

    expected_centre = ...
        x0 + velocity_values(u_index)* ...
        snapshot_times(validation_time_index);

    xlim([ ...
        max(0,expected_centre-90), ...
        min(channel_length,expected_centre+90)]);

    xlabel('x (m)');
    ylabel('C (kg/m)');
    title(sprintf('u = %.4f m/s',velocity_values(u_index)));

    grid on;

    if u_index == 1
        legend('LBM Numerical','Analytical');
    end

    hold off;
end

sgtitle('Numerical vs Analytical at t = 20000 s');

exportgraphics( ...
    gcf, ...
    fullfile(results_directory, ...
    '1D_velocity_validation_20000s.png'), ...
    'Resolution',300);

% Centroid position at 20000 s
selected_time = 20000;
selected_rows = Time_s == selected_time;

figure;

plot( ...
    Velocity(selected_rows), ...
    NumericalCentroidX_m(selected_rows), ...
    'o-','LineWidth',1.5);

hold on;

plot( ...
    Velocity(selected_rows), ...
    ExpectedCentroidX_m(selected_rows), ...
    '--','LineWidth',1.5);

xlabel('Velocity u (m/s)');
ylabel('Plume centroid x (m)');
title('Effect of Velocity on Plume Position at t = 20000 s');

legend('LBM Numerical','Theoretical');

grid on;
hold off;

exportgraphics( ...
    gcf, ...
    fullfile(results_directory, ...
    '1D_velocity_centroid_20000s.png'), ...
    'Resolution',300);

fprintf('\nVelocity experiment completed.\n');
fprintf('Results saved in: %s\n', results_directory);