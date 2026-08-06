%% Plot comparison - Proposed Case

% Simulation Case
sim = 3; % 1.REF-PCM / 2.REF-CPCM / 3.REF-DOUBLE

% Save figure plots
save = false;

% Results folders
fluent_folder = 'C:\Users\renat\Desktop\Fluent\Reference\Double\DoublePouch_files\dp0\FFF\Fluent';
vcm_folder = 'C:\Users\renat\Documents\SVN\VCM-PCM\trunk\Results\Results_trunk_DOUBLE_Mesh4_PCM';
save_folder = 'C:\Users\renat\Desktop\Fotos';

%% Reference Files
[dt, tmax, xx, time, prefix] = select_reference(sim);

%% Read Fluent data
addpath(genpath(fluent_folder));
max_temp_file = which('\max-batt-temp-rfile.out');
min_temp_file = which('\min-batt-temp-rfile.out');
liq_frac_file = which('\ave-liq-frac-rfile.out');

if ~exist(max_temp_file, 'file') || ~exist(min_temp_file, 'file') 
    max_temp_out = [fluent_folder '\max-temp-rfile.out'];
    min_temp_out = [fluent_folder '\min-temp-rfile.out'];
    copyfile(max_temp_out, max_temp_file);
    copyfile(min_temp_out, min_temp_file);
end

% Max Temp
opts = detectImportOptions(max_temp_file, 'FileType','text');
opts.VariableNames = ["Time", "Tmax", "FlowTime"];
data_max = readtable(max_temp_file, opts);
% Min Temp
opts = detectImportOptions(min_temp_file, 'FileType','text');
opts.VariableNames = ["Time", "Tmin", "FlowTime"];
data_min = readtable(min_temp_file, opts);
% Liq Fraction
opts = detectImportOptions(liq_frac_file, 'FileType','text');
opts.VariableNames = ["Time", "F", "FlowTime"];
data_f = readtable(liq_frac_file, opts);

rmpath(genpath(fluent_folder));

%% Read code data
addpath(genpath(vcm_folder));
log_file = which('\Log.dat');
opts = detectImportOptions(log_file);
opts.VariableNames = ["Time", "Tmax","Tmin", "DeltaT", "F"];
data = readtable(log_file, opts);

rmpath(genpath(vcm_folder));

close all
fig = figure(1);
plot(data.Time, data.Tmin, '-')
hold on
plot(data_min.FlowTime, data_min.Tmin - 273.15, '-')
hold on
yy  = spline(tmax.time, tmax.tmax, xx) - spline(dt.time, dt.dt, xx);
plot(xx, yy, 's')
grid on
%title('Minimum Temperature in Battery Module')
xlabel('Time (s)')
ylabel('Temperature (°C)')
%legend('VCM', 'Fluent', 'Reference', 'Location', 'northwest')
legend('VCM', 'Fluent', 'Reference', 'Location', 'northwest')
axis([0 500 25 45]);
if save
    exportgraphics(fig, [save_folder '\' prefix 'Comparison_MinTemp.pdf'], 'ContentType', 'vector');
end

fig = figure(2);
plot(data.Time, data.Tmax, '-')
% hold on
% plot(data_max.FlowTime, data_max.Tmax - 273.15, '-')
hold on
yy  = spline(tmax.time, tmax.tmax, xx);
plot(xx, yy, 's')
grid on
%title('Maximum Temperature in Battery Module')
xlabel('Time (s)')
ylabel('Temperature (°C)')
%legend('VCM', 'Fluent', 'Reference', 'Location', 'northwest')
legend('VCM', 'Reference', 'Location', 'northwest')
axis([0 500 25 45]);
if save
    exportgraphics(fig, [save_folder '\' prefix 'Comparison_MaxTemp.pdf'], 'ContentType', 'vector');
end

filter_data = (mod(data_f.FlowTime, 10) == 0);
fig = figure(3);
plot(data.Time, data.F, '-')
hold on
plot(data_f.FlowTime(filter_data), data_f.F(filter_data), 's')
grid on
%title('Average Liquid Fraction in Battery Module')
ylabel('Liquid Fraction')
%legend('VCM', 'Fluent', 'Reference', 'Location', 'northwest')
legend('VCM', 'Fluent', 'Location', 'northwest')
axis([0 500 0 1]);
if save
    exportgraphics(fig, [save_folder '\' prefix 'Comparison_LiqFrac.pdf'],  'ContentType', 'vector');
end

%% RMSE TMax
xx = 0:0.001:500.001;
yy = spline(tmax.time, tmax.tmax, xx);

soma = 0;
diff_max = 0;
for i = 1:numel(data.Tmax)
    diff = (data.Tmax(i) - yy(i))^2;
    soma = diff + soma;
    if abs(data.Tmax(i) - yy(i)) > diff_max
        diff_max =  abs(data.Tmax(i) - yy(i));
    end
end
rmse_tmax = soma/i;
diff_tmax = diff_max;

%% RMSE TMin
yy  = spline(tmax.time, tmax.tmax, xx) - spline(dt.time, dt.dt, xx);

soma = 0;
diff_max = 0;
for i = 1:numel(data.Tmin)
    diff = (data.Tmin(i) - yy(i))^2;
    soma = diff + soma;
    if abs(data.Tmin(i) - yy(i)) > diff_max
        diff_max = abs(data.Tmin(i) - yy(i));
    end
end
rmse_tmin = soma/i;
diff_tmin = diff_max;

%% RMSE LiqFrac
xx = 0:0.001:500.001;
yy  = spline(data_f.FlowTime, data_f.F, xx);

soma = 0;
diff_max = 0;
for i = 1:numel(data.F)
    diff = (data.F(i) - yy(i))^2;
    soma = diff + soma;
    if abs(data.F(i) - yy(i)) > diff_max
        diff_max = abs(data.F(i) - yy(i));
    end
end
rmse_f = soma/i;
diff_f = diff_max;
