%% Plot comparison - Proposed Case

% Simulation Case
sim = 3; % 1.REF-PCM / 2.REF-CPCM / 3.REF-DOUBLE / 4.3C-PCM / 5.3C-CPCM / 6.3C-DOUBLE

% Save figure plots
save = true;

% Results folders
fluent_folder = 'Fluent_Double_h50_dx_1p0_Mushy';
vcm_folder = 'Results_2024-11-17_10-28_REV21_1p0_resmax_10-3';
    
%% Read Fluent data
[filepath,~,~] = fileparts(mfilename('fullpath'));
results_fullpath = [extractBefore(filepath, '\Validation') '\Results\' fluent_folder];
max_temp_file = [results_fullpath '\max-temp-rfile.txt'];
min_temp_file = [results_fullpath '\min-temp-rfile.txt'];

if ~exist(max_temp_file, 'file') || ~exist(min_temp_file, 'file') 
    max_temp_out = [results_fullpath '\max-temp-rfile.out'];
    min_temp_out = [results_fullpath '\min-temp-rfile.out'];
    copyfile(max_temp_out, max_temp_file);
    copyfile(min_temp_out, min_temp_file);
end

% Max Temp
opts = detectImportOptions(max_temp_file);
opts.VariableNames = ["Time", "Tmax", "FlowTime"];
data_max = readtable(max_temp_file, opts);
% Min Temp
opts = detectImportOptions(min_temp_file);
opts.VariableNames = ["Time", "Tmin", "FlowTime"];
data_min = readtable(min_temp_file, opts);

%% Read code data 
[filepath,~,~] = fileparts(mfilename('fullpath'));
results_fullpath = [extractBefore(filepath, '\Validation') '\Results\' vcm_folder];
opts = detectImportOptions([results_fullpath '\Log.dat']);
opts.VariableNames = ["Time", "Tmax","Tmin", "DeltaT", "F"];
data = readtable([results_fullpath '\Log.dat'], opts);

%close all
fig = figure(1);
plot(data.Time, data.Tmin, '-')
hold on
plot(data_min.Time, data_min.Tmin - 273.15, '-')
grid on
title('Minimum Temperature in Battery Module')
xlabel('Time (s)')
ylabel('Temperature (°C)')
legend('Code','Fluent', 'Location', 'northwest')
axis([0 500 25 45]);
if save
    saveas(fig, [fluent_folder '_MinTemp.png'])
end

fig = figure(2);
plot(data.Time, data.Tmax, '-')
hold on
plot(data_max.Time, data_max.Tmax - 273.15, '-')
grid on
title('Maximum Temperature in Battery Module')
xlabel('Time (s)')
ylabel('Temperature (°C)')
legend('Max Temp - Code','Max Temp - Fluent')
axis([0 500 25 45]);
if save
    saveas(fig, [fluent_folder '_MaxTemp.png'])
end

fig = figure(3);
plot(data.Time, data.Tmax - data.Tmin, '-')
hold on
plot(data_min.Time, data_max.Tmax - data_min.Tmin, '-')
grid on
title('Maximum Temperature Difference in Battery Module')
ylabel('Temperature difference (°C)')
legend('\DeltaTemp - Code','\DeltaTemp - Fluent', 'Location', 'northwest')
axis([0 500 0 10]);
if save
    saveas(fig, [fluent_folder '_DeltaTemp.png'])
end