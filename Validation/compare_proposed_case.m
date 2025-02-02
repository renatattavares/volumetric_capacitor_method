%% Plot comparison - Proposed Case

% Simulation Case
sim = 3; % 1.REF-PCM / 2.REF-CPCM / 3.REF-DOUBLE / 4.3C-PCM / 5.3C-CPCM / 6.3C-DOUBLE

% Save figure plots
save = true;

% Results folder
results_folder = 'Results_2024_12_15_15_07_14';
    
%% Result files
[filepath,~,~] = fileparts(mfilename('fullpath'));
results_fullpath = [extractBefore(filepath, '\Validation') '\Results\' results_folder];
opts = detectImportOptions([results_fullpath '\Log.dat']);
opts.VariableNames = ["Time", "Tmax","Tmin", "DeltaT", "F"];
data = readtable([results_fullpath '\Log.dat'], opts);
  
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

%{
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

%close all
fig = figure(3);
yy  = spline(tmax.time, tmax.tmax, xx) - spline(dt.time, dt.dt, xx);
plot(xx, yy, 's')
hold on
plot(data.Time, data.Tmin, '-')
grid on
title('Minimum Temperature in Battery Module')
xlabel('Time (s)')
ylabel('Temperature (°C)')
legend('Reference','Result', 'Location', 'northwest')
axis([0 time 25 45]);
if save
    saveas(fig, [results_folder '_MinTemp.png'])
end
fig = figure(4);
ax = axes;
yyaxis left
yy  = spline(tmax.time, tmax.tmax, xx);
plot(xx, yy, 's')
hold on
plot(data.Time, data.Tmax, '-')
grid on
xlabel('Time (s)')
ylabel('Temperature (°C)')
axis([0 time 25 45]);

yyaxis right
yy = spline(dt.time, dt.dt, xx);
plot(xx, yy, 's')
hold on
plot(data.Time, data.DeltaT, '-')
grid on
ylabel('Temperature difference (°C)')
legend('Max Temp - Reference','Max Temp - Result', '\DeltaTemp - Reference','\DeltaTemp - Result', 'Location', 'northwest')
axis([0 time 0 10]);
if save
    saveas(fig, [results_folder '_MaxTemp_DeltaTemp.png'])
end

% Set the color of each axis to black
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

%% RMSE
try
    ref = spline(tmax_ref.time, tmax_ref.Tmax, xx);
    mine = spline(data.Time, data.Tmax, xx);
    err = sqrt(sum((ref(:)-mine(:)).^2) / numel(ref));
catch ME
   %warning(ME.message) 
end

%ref = spline(tmax_ref.time, tmax_ref.Tmax, xx) - spline(deltat_ref.time, deltat_ref.DeltaT, xx);
%mine = spline(data.Time, data.Tmin, xx);
%err = sqrt(sum((ref(:)-mine(:)).^2) / numel(ref));

%}