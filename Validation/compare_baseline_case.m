%% Plot comparison - Baseline Case

% Results folder
%filename = 'C:\Users\renat\Documents\SVN\VCM-PCM\branches\2024_06_07_ModelValidation_PouchCells\Results\Baseline_0p0005';
%filename = 'C:\Users\renat\Documents\SVN\VCM-PCM\trunk\Results\Results_2024-08-04_00-57';
filename = 'C:\Users\renat\Documents\SVN\VCM-PCM\branches\2024_08_15_NewNavigation\Results\Results_2024-08-18_16-48';

% Reference Files
deltat_ref = readtable('DeltaT_ref.txt');
tmax_ref = readtable('Tmax_ref.txt');

% Result files
opts = detectImportOptions([filename '\Log.dat']);
opts.VariableNames = ["Time", "Tmax","Tmin", "DeltaT", "F"];
data = readtable([filename '\Log.dat'], opts);

% Print titles
titles = 1;

%% Create plots
close all
figure(1)
xx = 0:50:1200;
yy  = spline(tmax_ref.time, tmax_ref.Tmax, xx);
plot(xx, yy, 'o')
hold on
plot(data.Time, data.Tmax)
grid on
xlabel('Time (s)')
ylabel('Temperature (°C)')
axis([0 1300 25 60]);
legend('Reference Case','Simulation Result', 'Location', 'northwest')
if titles
    title('Maximum Temperature in Battery Module')
end

figure(2)
xx = 0:50:1200;
yy  = spline(deltat_ref.time, deltat_ref.deltat, xx);
plot(xx, yy, 'o')
hold on
plot(data.Time, data.DeltaT)
grid on

xlabel('Time (s)')
ylabel('Temperature difference (°C)')
axis([0 1300 0 20]);
legend('Reference Case','Simulation Result', 'Location', 'northwest')
if titles
    title('Maximum Temperature Difference in Battery Module')
end

figure(3)
xx = 0:60:1200;
yy  = spline(tmax_ref.time, tmax_ref.Tmax, xx) - spline(deltat_ref.time, deltat_ref.deltat, xx);
plot(xx, yy, 'o')
hold on
plot(data.Time, data.Tmin)
grid on
xlabel('Time (s)')
ylabel('Temperature (°C)')
axis([0 1300 25 50]);
legend('Reference','Simulation Result', 'Location', 'northwest')
if titles
    title('Minimum Temperature in Battery Module');
end

figure(4)
ax = axes;
yyaxis left
xx = 0:50:1200;
yy  = spline(tmax_ref.time, tmax_ref.Tmax, xx);
plot(xx, yy, 'o')
hold on
plot(data.Time, data.Tmax, '-')
grid on
xlabel('Time (s)')
ylabel('Temperature (°C)')
axis([0 1300 25 60]);

yyaxis right
yy = spline(deltat_ref.time, deltat_ref.deltat, xx);
plot(xx, yy, 'o')
hold on
plot(data.Time, data.DeltaT, '-')
grid on
ylabel('Temperature difference (°C)')
axis([0 1300 0 30]);
legend('Maximum Temperature - Reference','Maximum Temperature - Simulation Result', '\DeltaTemperature - Reference','\DeltaTemperature - Simulation Result', 'Location', 'northwest')

% Set the color of each axis to black
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

%% RMSE
xx = 0:50:1200;
ref = spline(tmax_ref.time, tmax_ref.Tmax, xx);
mine = spline(data.Time, data.Tmax, xx);
err = sqrt(sum((ref(:)-mine(:)).^2) / numel(ref));

