%  Input info
PlotsPath = uigetdir('C:\Users\renat\Documents\SVN\VCM-PCM\Results');
prompt = {'Enter simulation final time:'};
dlgtitle = 'Input';
fieldsize = [1 45];
definput = {'300'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
FinalTime = str2num(answer{1});
Tinitial = 20;
Nx = 200;
Ny = 200;
Lx = 5; % [cm]
Ly = 5; % [cm]
dx = Lx/Nx; 
dy = Ly/Ny;
mesh = zeros(Nx*Ny, 2);

% Build mesh reference
id = 1;
for i = 1:Nx
   
    x = (i-1)*dx + 0.5 * dx;
    
    for j = 1:Ny
        y = (j-1)*dy + 0.5 * dy;
        mesh(id,1) = x;
        mesh(id,2) = y;
        id = id+1;
    end
    
end
    
% Filter files
cd(PlotsPath)
DirStruct = dir(PlotsPath);
DirCell = struct2cell(DirStruct);
DirFiles = DirCell(1,:)';
files = contains(DirFiles, ".vtk");
files = DirFiles(files);
prop = extractBefore(files{1}, "_");

% Get data range
max_value = 35;
if strcmp(prop, 'Temp')
    min_value = Tinitial;
    unit = "[°C]";
else
    min_value = 0;
    unit = "";
end

% Create Pics folder
pics_folder = PlotsPath + "\Plots";
if ~exist(pics_folder, 'dir')
    try
        mkdir(pics_folder);
    catch ME
        rethrow(ME);
    end
end

% Generate pics
for i = 1:length(files)
   
    filename = files{i};
    time = extractAfter(extractBefore(filename, '.'), '_');
    data = fileread(filename);
    data = split(extractAfter(data, 'LOOKUP_TABLE my_table'));
    data = data(~cellfun(@isempty, data));
    data = cellfun(@str2num, data);
    
    % Plot temperature
    s = scatter(mesh(:,1), mesh(:,2), 'o', 'filled');
    s.CData = data;
    set(gcf, 'Position', [500,200,500,400]);
    figName = prop + "in " + time + " secs" + unit;
    title(figName)
    xlabel('X Coordinates [cm]')
    ylabel('Y Coordinates [cm]')
    c = colorbar;
    if strcmp(prop, 'Temp')
        %c.Limits = [Tinitial ceil(max_value)];
        %c.Ticks = linspace(Tinitial, ceil(max_value), 10);
        %c.TickLabels = round(linspace(Tinitial, ceil(max_value), 10));
    else
        %c.Limits = [0 1];
    end
    saveas(gcf, pics_folder + '\' + prop + '_' + time + ".png")
    close all

end