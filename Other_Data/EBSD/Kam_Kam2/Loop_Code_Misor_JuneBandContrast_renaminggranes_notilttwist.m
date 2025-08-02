clear all 
close all

%% ADD MTEC ptah
addpath("C:\Users\charl\Documents\Programs\mtex-5.11.2");
startup_mtex;

%% Specify Crystal and Specimen Symmetries _ Load depending on analysis setup
% crystal symmetry
% This is used to ID subgrains
subgrain_angle=0.3; 
% This is used to pull out what counts as a separate grain when plotting each grain separatly
grain_div_angle=1; 
color_angle=3;
len_filter=1 % filter out shorter grain boundaries than this um
results = {};

CS = {... 
  'notIndexed',...
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [8.358 8.358 8.358], 'mineral', 'Chromite', 'color', 'light green'),...
  crystalSymmetry('Pbca', [7.048 7.193 18.159], 'mineral', 'Enstatite', 'color', 'light yellow')};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
%% Loading individual files

day_set='CD'
% path to filesf
if strcmp(day_set, 'comp');
    pname = 'C:\Users\penny\Box\Berkeley_new\EBSD\Detector_Comparison'
elseif strcmp(day_set, 'day3_ctf');
    pname = 'C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Day three EBSD\EBSD ML Day 3\ctf_files_day3';
elseif strcmp(day_set, 'day3');
    pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Day three EBSD\EBSD ML Day 3\h5oina';
elseif strcmp(day_set, 'day2');
    pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Day2_MLB\Data_for_paper_notests';
elseif strcmp(day_set, 'day1');
    pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\OldEBSDProjects_beforeReexport\Data_ML_EBSD\Day1';
elseif strcmp(day_set, 'day4');
       pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\day four ebsd\Day 4 ebsd\ctf';
elseif strcmp(day_set, 'day5');
       pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Day five Christina + ML\Ori1';
       CS = {... 
  'notIndexed',...
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [8.4 8.4 8.4], 'mineral', 'Chromite', 'color', [0.56 0.74 0.56]),...
  crystalSymmetry('-3m1', [4.9 4.9 5.4], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Quartz', 'color', [0.85 0.65 0.13]),...
  crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite  Opx AV77', 'color', [0.94 0.5 0.5]),...
  crystalSymmetry('12/m1', [9.7 9 5.3], [90,105.63,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Diopside   CaMgSi2O6', 'color', [0 0 0.55])};



elseif strcmp(day_set, 'CD');
    pname="P:\WORK-GENERAL\POSTDOC-UCB\BERKELEY-VIBE\Documents\Projects\Kamaehu2024\2large4GIT\EBSD\Kam_Kam2\Kam_Kam2\ctf"
    CS = {... 
  'notIndexed',...
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [8.4 8.4 8.4], 'mineral', 'Chromite', 'color', [0.56 0.74 0.56]),...
  crystalSymmetry('-3m1', [4.9 4.9 5.4], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Quartz', 'color', [0.85 0.65 0.13]),...
  crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite  Opx AV77', 'color', [0.94 0.5 0.5]),...
  crystalSymmetry('12/m1', [9.7 9 5.3], [90,105.63,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Diopside   CaMgSi2O6', 'color', [0 0 0.55])};

elseif strcmp(day_set, 'CC');
     pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Christina_Code'

CS = {... 
  'notIndexed',...
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [8.4 8.4 8.4], 'mineral', 'Chromite', 'color', [0.56 0.74 0.56])};
       %pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\day four ebsd\Day 4 ebsd\h5oina';
elseif strcmp(day_set, 'CC_Penn');
pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Christina_Code\PennState'
     CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [8.4 8.4 8.4], 'mineral', 'Chromite', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.56 0.74 0.56])};

elseif strcmp(day_set, 'day6');
    pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Day6_EBSD\ctffiles'

elseif strcmp(day_set, 'Kil')
    pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Kilauea_EBSD_for_Comparison'

else
    error('Invalid day_set value');
end

folder = fullfile(pname, 'mis2mean_unfilt');

% Load Individual EBSD Map  - special cases, MLP_9, MP2_63, MLP_43,
% MLP_39,MLP_49, MP2_67, M, MP2_58B_hough60_refined 
% Get a list of all .ctf files in the directory
%files = dir(fullfile(pname, '\MLP_9 - EBSD Data.ctf')); % was '*.ctf'
files = dir(fullfile(pname, '\*.ctf')); % was '*.ctf'#
%files = dir(fullfile(pname, '\Day 4 ebsd ori2 ori2_ol1 Map Data 15.h5oina')); % was '*.ctf'#

for i = 1:length(files)
    display(files(i))
% We need to clear variables each time thro





if contains(files(i).name, 'ctf')
    filename = ['\' files(i).name];
    filename=strrep(filename, ' - EBSD Data','');
    filename=strrep(filename, '.ctf','');
    fname = fullfile(pname, files(i).name);

    % Taking from Charlie - matches everything in crystal 
    %ebsd = EBSD.load(fname,CS,'interface','ctf','convertEuler2SpatialReferenceFrame');
    ebsd = EBSD.load(fname, 'convertEuler2SpatialReferenceFrame', 'interface', 'ctf');
    % rotate everything 180 degrees around x to get the correct map
    ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
    % then rotate only the orientations 180 degrees around z
    ebsd = rotate(ebsd,rotation('axis',zvector,'angle',180*degree),'keepXY');
%ebsd=rotate(ebsd, 180*degree, 'keepXY')



else
    
    filename = ['\' files(i).name];
    filename=strrep(filename, ' Map Data ','_');
    filename=strrep(filename, '.h5oina','');
    fname = fullfile(pname, files(i).name);
% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','h5oina', 'convertEuler2SpatialReferenceFrame');
% rotate everything 180 degrees around x to get the correct map
ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
ebsd = rotate(ebsd,rotation('axis',yvector,'angle',180*degree),'keepXY');

end
ebsd.unitCell = ebsd.unitCell;


%% Data Filtering - remove misindexed pixels
[grains, ebsd.grainId, ebsd.mis2mean]= calcGrains(ebsd,'angle',10*degree);
ebsd(grains(grains.grainSize<2)).phase=0; %removes misindexed pixels 


%% Apply 5 neighbour Kuwahara filter and fill grains
F = KuwaharaFilter; 
F.numNeighbours = 5;
ebsd= smooth(ebsd('indexed'), F, 'fill',grains);




%% Subgrain reconstruction, puts true grain boundaries at grain_div_angle, subgrain boundaries at subgrain_angle
[grains,ebsd.grainId,ebsd.mis2mean]= calcGrains(ebsd,'angle',[subgrain_angle*degree,grain_div_angle*degree]); % Want this big, so it pulls out each grain separatly. 

% loops over smoothing multiple times.
iter=[1,5,8];
for Z=1:length(iter)
    grains=smooth(grains, iter(Z));
end
grains=grains(grains.area >5000); % grains(grains.grainSize > max(grains.grainSize/30))  % grains(grains.area > 100); % Only use grains >100 um^2 in size (so 10X 10 um) % 
grains=grains('Fo'); %Construct grains variable that is restricted to only the phase forsterite

remainingGrainIds = [grains.id];
isRemainingGrain = ismember(ebsd.grainId, remainingGrainIds);
ebsd_Fo = ebsd(isRemainingGrain);



%% Mis 2 mean plot
mtexFig = newMtexFigure('layout',[1, 3]); % 1 row, 2 columns
% First axis - showing band contrast
ax1 = mtexFig.nextAxis;
plot(ebsd_Fo, ebsd_Fo.bc);
colormap(ax1, gray);
title(ax1, 'Band contrast');
hold(ax1, 'on');

% Second axis
ax2 = mtexFig.nextAxis;
plot(ebsd_Fo, ebsd_Fo.mis2mean.angle ./ degree);
colormap(ax2, jet);
title(ax2, 'Mis2Mean Coloring');
hold(ax2, 'on');

% Plot additional elements on the second axis
plot(grains.innerBoundary, 'lineColor', 'k', 'width', 2); %, 'DisplayName', 'subgrain inner'
plot(grains.boundary, 'lineColor', 'r'); % 'DisplayName', 'Grain boundary', 'fontsize', 5

% Set color limits for the second axis
clim(ax2, [0 color_angle]);

% Create and configure the colorbar for the second axis
h = colorbar('peer', ax2); % Use MATLAB's colorbar function
h.Label.String = 'Mis2Mean angle';  % Set the label text


% Loop through grains to label them with their GOS value. 
for j = 1:length(grains)
    
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grains(j).boundary.x);
    centerY = mean(grains(j).boundary.y);
    
    % Get the grain ID and GOS value
    grainId = grains(j).id;
    GOS = grains(j).GOS ./ degree;
    
    
    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('ID: %d\nGOS: %.2f', grainId, GOS), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 1);
end




% Lets add the IPF figure --------------
mtexFig.nextAxis;
title('IPF Coloring')

% Lopp through grains to plot an IPF figure to show misorientation angle +
% direction relative to mean grain 
for z = 1:length(grains)
    
    % Get EBSD data for the current grain
    currentGrainIdsz = grains(z).id;  % IDs of pixels in the current grain

     % Combine grain ID with filename
    namez = [num2str(currentGrainIdsz) '_' filename];
    currentGrainEBSDz = ebsd_Fo(ebsd_Fo.grainId == currentGrainIdsz);  % Subset EBSD data for current grain
    oM = ipfHSVKey(currentGrainEBSDz);
    oM.inversePoleFigureDirection = orientation(grains(z).meanOrientation) * oM.whiteCenter;
    oM.colorPostRotation = reflection(yvector);
    oM.maxAngle = 2*degree;
    color=oM.orientation2color(currentGrainEBSDz.orientations);
    hold on
    [~,mP]=plot(currentGrainEBSDz, color);


    
    
end

% Check if the folder 'mis2mean_unfilt' exists in the specified path; if not, create it

if ~exist(folder, 'dir')
    mkdir(folder);
end

savePath = fullfile(folder, [strrep(filename, '.ctf', ''), '_mis2mean.png']);
saveas(gcf, savePath)

close all 

%-------------------------------------------------------------------------
% Start afresh with a new figure! This is to show each grain, and its
% misorientation signature 
%------------------------------------------------------------------------
mtexFig = newMtexFigure('layout',[2, length(grains)]); % 2 row, N columns
figure(mtexFig.parent);  % Ensure we are modifying the correct figure
set(mtexFig.parent, 'Position', [100, 100, 1200, 1200]);  % Adjust these values as needed
hold on


% Loop through each grain again for IPF coloring. 
for k = 1:length(grains)
    
    x=1; % x value for plotting
    y=k; % column based on grain
    
    % Get EBSD data for the current grain
    currentGrainIds = grains(k).id;  % IDs of pixels in the current grain
    display(currentGrainIds)
     % Combine grain ID with filename
    name = [num2str(currentGrainIds) '_' filename];
    GOS = grains(k).GOS ./ degree;



    currentGrainEBSD = ebsd_Fo(ebsd_Fo.grainId == currentGrainIds);  % Subset EBSD data for current grain
    oM = ipfHSVKey(currentGrainEBSD);
    oM.inversePoleFigureDirection = grains(k).meanOrientation * oM.whiteCenter; % Previously was orientation(grains(k).meanOrientation) * oM.whiteCenter;
    oM.maxAngle = grain_div_angle*degree;
    color=oM.orientation2color(currentGrainEBSD.orientations);

    % Lets add the grain number
    
    

    % This is plotting the pretty IPF map. 
    mtexFig.nextAxis(x, y)
    hold on
    [~,mP]=plot(currentGrainEBSD, color);
    centerX = mean(grains(k).boundary.x);
    centerY = mean(grains(k).boundary.y);
    text(centerX, centerY, sprintf('ID: %d\nGOS: %.2f', currentGrainIds , GOS), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 1);
   


    % Lets calculate grain reference orientation deviation (GROD) 
    % the misorientation between each pixel orientation to the grain mean orientation
    ori=currentGrainEBSD.orientations;
    mis2mean=inv(grains(k).meanOrientation) .* ori;
    GOS=currentGrainEBSD.grainMean(mis2mean.angle, currentGrainIds);
    


     gb_C_m=[grains(k).innerBoundary];

    if isempty(gb_C_m)
        disp('Grain boundary is empty!');
         x=2; %update to move onto the next figure 
        mtexFig.nextAxis(x, y);
        results = [results; {filename, grains(k).id, grains(k).GOS./degree, grains(k).grainSize, ...
    grains(k).area, grains(k).subBoundaryLength, grains(k).equivalentRadius, 0, 0,0, 0, 0, 0, 0,subgrain_angle, grain_div_angle  }];

        continue; % Skip the current iteration if there is no boundary data
    end


   

%----------------------- First, lets filter boundaries based on
%misorientation angles ------------------------------------------%
         
% 
    uppangle = subgrain_angle; % Plots misorientation axes that have angles greater than this
    ang = 10; % Plots misorientation axes that have angles less than this
% 
%     % Conditions
    cond1 = gb_C_m.misorientation.angle./degree > uppangle;
    cond2 = gb_C_m.misorientation.angle./degree < ang;
    cond3 = gb_C_m.misorientation.angle./degree > ang;

    cond=uppangle < gb_C_m.misorientation.angle./degree & gb_C_m.misorientation.angle./degree < ang;

    condlength=gb_C_m.segLength>len_filter;
    

        if ~any(cond&condlength)
        disp('No boundaries wit selected angles');
        x=2; %update to move onto the next figure 
        mtexFig.nextAxis(x, y);
        continue; % Skip the current iteration if there is no boundary data
        end


    plot(gb_C_m(condlength&cond), 'linecolor','k','linewidth',1.5) % This is all boundaries. 
    %plot(grains.innerBoundary, 'lineColor', 'y', label='subgrain inner')
    mP.micronBar.visible = 'off';

    % Title for each subplot
    formattedFilename = strrep(filename, '\', '\\');
formattedFilename = strrep(formattedFilename, '_', '-');

% Set the title with specified font size
title(['Grain ', num2str(grains(k).id), ' from ', formattedFilename], 'FontSize', 10);






Tot_length=sum(gb_C_m.segLength);
AvRotAxis=mean(gb_C_m.misorientation(cond&condlength));
prop_length=Tot_length/((grains(k).grainSize)^0.5);

results = [results; {filename, grains(k).id, grains(k).GOS./degree, grains(k).grainSize, ...
    grains(k).area, grains(k).subBoundaryLength, grains(k).equivalentRadius, prop_length, Tot_length,...
    0, 0, 0, 0, 0, subgrain_angle, grain_div_angle  }];


    % Plot data on next axis
    x=2; %update to move onto the next figure 
    mtexFig.nextAxis(x, y); % Now move onto the next axis - to see what you get IPF wise 
    hold on
    plotAxisDistribution(gb_C_m.misorientation(cond&condlength), 'contourf');

    % Lets save these axes
    folder_mis = fullfile(pname, 'Mis_axes');
    if ~exist(folder_mis, 'dir')
        mkdir(folder_mis);
    end

    %if (grains(k).GOS./degree >0.25)
    
    savename=sprintf('%s_grain%d_misorientation', filename, grains(k).id);
    save_mis=gb_C_m.misorientation;
    fullFilePath_axis=fullfile(folder_mis, savename);
    save(fullFilePath_axis, 'save_mis')
    %end
    %save('')




end
% Save the entire figure with all subplots
savePath = fullfile(folder, [filename, '_IPF_All_Grains.png']);
mkdir(folder); % Ensure the folder exists
saveas(gcf, savePath);




close all

varsToKeep = {'files', 'folder', 'pname', 'len_filter', 'CS', 'i', 'subgrain_angle', 'grain_div_angle', 'color_angle', 'results'};
allVars = who;% Get a list of all variables in the workspace


varsToClear = setdiff(allVars, varsToKeep); % Determine which variables to clear

clear(varsToClear{:}); % Clear the selected variables


end
%% 

%%
fullmatrix=fullfile(pname, 'test.xlsx');
headers={'Filename', 'grainID',	'GOS' 'Grain Size (pixels)', 'Grain Size (um2)',...
    'subBoundaryLength (um)', 'equivalentRadius (um)',	'GB Length/Sqrt Size', 'Total GB length',	...
    'Tilt length',	'Twist length',	'Perc Tilt',	'Perc Twist', 'Perc unclassified', 'subgrain_angle', 'grain_div_angle'};
% results = [results; {filename, grains(k).id, grains(k).GOS./degree, grains(k).grainSize, ...
%     grains(k).area, grains(k).subBoundaryLength, grains(k).equivalentRadius, prop_length, Tot_length,0, 0, 0, 0, 0, subgrain_angle, grain_div_angle  }];

results_with_head=[headers; results];
writecell(results_with_head, fullmatrix)

% %% Load one as a test
% % Set the folder path
% folder_path = fullfile(folder, 'Mis_axes');
% 
% % List all files in the folder
% file_list = dir(folder_path);
% 
% % Initialize an empty cell array to store the loaded data
% data_combined = [];
% 
% % Loop through each file in the folder
% for i = 1:length(file_list)
%     % Get the current file name
%     current_file = file_list(i).name;
% 
%     % Check if the current file contains 'MLP_43' in its name
%     if contains(current_file, 'MLP_43') && strcmp(current_file(end-3:end), '.mat')
%         % Construct the full path to the .mat file
%         file_path = fullfile(folder_path, current_file);
% 
%         % Load the specific variable from the file
%         data = load(file_path, 'save_mis');
% 
%         % Append the loaded data to the combined array
%         data_combined = [data_combined; data.save_mis];
%     end
% end
% 
% % Plot the combined data
% plotAxisDistribution(data_combined, 'contourf');


%% Now Load them all


% % Set the folder path
% folder_path = [folder '/Mis_axes'];
% 
% % Get a list of all .mat files in the folder
% file_list = dir(fullfile(folder_path, '*.mat'));
% 
% % Initialize a cell array to store the loaded data
% data_combined = [];
% 
% % Loop through each file and load the data
% for i = 1:length(file_list)
%     file_path = fullfile(folder_path, file_list(i).name);
%     data = load(file_path, 'save_mis'); % Load the specific variable from each file
%     data_combined = [data_combined; data.save_mis]; % Append to the combined variable
% end
% 
% % Plot the combined data
% plotAxisDistribution(data_combined, 'contourf');


  