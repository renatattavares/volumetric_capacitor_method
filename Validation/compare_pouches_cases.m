%% Plot comparison - Proposed Case

% Save figure plots
save = false;

% Results folder
PCM_results_folder = 'Results_2024-11-17_12-12_REV21_FullPCM_dx_0p5';
CPCM_results_folder = 'Results_2024-11-20_10-57_REV21_FullCPCM_dx_1p0';
DOUBLE_results_folder = 'Results_2024-11-16_11-41_REV21_0p5_resmax_10-4';

%% Result files
[filepath,~,~] = fileparts(mfilename('fullpath'));
PCM_results_fullpath = [extractBefore(filepath, '\Validation') '\Results\' PCM_results_folder];
CPCM_results_fullpath = [extractBefore(filepath, '\Validation') '\Results\' CPCM_results_folder];
DOUBLE_results_fullpath = [extractBefore(filepath, '\Validation') '\Results\' DOUBLE_results_folder];
opts = detectImportOptions([PCM_results_fullpath '\Log.dat']);
opts.VariableNames = ["Time", "Tmax","Tmin", "DeltaT", "F"];

PCM_data = readtable([PCM_results_fullpath '\Log.dat'], opts);
CPCM_data = readtable([CPCM_results_fullpath '\Log.dat'], opts);
DOUBLE_data = readtable([DOUBLE_results_fullpath '\Log.dat'], opts);

%{
% Reference Files
switch sim
    case 1 % Reference - Pure PCM Pouch
        dt = readtable('REF_PCM_DT.txt');
        tmax = readtable('REF_PCM_TMAX.txt');
        xx = 0:20:500;
        time = 500;
    case 2 % Reference - Pure CPCM Pouch
        dt = readtable('REF_CPCM_DT.txt');
        tmax = readtable('REF_CPCM_TMAX.txt');
        xx = 0:20:500;
        time = 500;
    case 3 % Reference - Double Pouch
        dt = readtable('REF_DOUBLE_DT.txt');
        tmax = readtable('REF_DOUBLE_TMAX.txt');
        xx = 0:20:500;
        time = 500;
    case 4
        return
    case 5
        return
    case 6
        deltat_ref = readtable('DeltaT_proposed_3C.txt');
        tmax_ref = readtable('Tmax_proposed_3C.txt');
        xx = 0:50:1200;
        time = 1250;
end


 Create plots
%close all
% figure(1)
% yy  = spline(tmax_ref.time, tmax_ref.Tmax, xx);
% plot(xx, yy, 'o')
% hold on
% plot(data.Time, data.Tmax)
% grid on
% title('Maximum Temperature in Battery Module')
% xlabel('Time (s)')
% ylabel('Temperature (°C)')
% legend('Reference Case','Simulation Result', 'Location', 'northwest')
% axis([0 Time 25 45]);
% 
% figure(2)
% yy  = spline(deltat_ref.time, deltat_ref.DeltaT, xx);
% plot(xx, yy, 'o')
% hold on
% plot(data.Time, data.DeltaT)
% grid on
% title('Maximum Temperature Difference in Battery Module')
% xlabel('Time (s)')
% ylabel('Temperature difference (°C)')
% legend('Reference Case','Simulation Result', 'Location', 'northwest')
% axis([0 Time 0 10]);
%}

fig = figure(1);
ax = axes;
yyaxis left
plot(PCM_data.Time, PCM_data.Tmax, 'Color', 'red', 'LineWidth', 1)
hold on
plot(CPCM_data.Time, CPCM_data.Tmax, '-', 'Color', 	"#EDB120", 'LineWidth', 1)
hold on
plot(DOUBLE_data.Time, DOUBLE_data.Tmax, '-', 'Color', "#77AC30", 'LineWidth', 1)
grid on
xlabel('Time (s)')
ylabel('Temperature (°C)')
axis([0 time 25 45]);

yyaxis right
yy = spline(dt.time, dt.dt, xx);
plot(PCM_data.Time, PCM_data.DeltaT, '--', 'Color', 'red', 'LineWidth', 1)
hold on
plot(CPCM_data.Time, CPCM_data.DeltaT, '--', 'Color', 	"#EDB120", 'LineWidth', 1)
hold on
plot(DOUBLE_data.Time, DOUBLE_data.DeltaT, '--', 'Color', "#77AC30", 'LineWidth', 1)
grid on
ylabel('Temperature difference (°C)')
legend('PCM - Max Temp','CPCM - Max Temp', 'DOUBLE - Max Temp', 'PCM - \DeltaTemp','CPCM - \DeltaTemp', 'DOUBLE - \DeltaTemp', 'Location', 'northwest')
axis([0 time 0 10]);
if save
    saveas(fig, ['AllPouches_MaxTemp_DeltaTemp.png'])
end

