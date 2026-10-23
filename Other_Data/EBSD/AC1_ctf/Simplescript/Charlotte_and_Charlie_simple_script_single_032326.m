%% Charlotte Devitre and Charlie Gordon's map script
% If you use this script, please give us both credit and refer to the github. 
% Updated June 2025 for mtex 6.0.0 (using updated grain calculation code)

close all
clear

% Start up mtex
addpath("/Users/cljd/Documents/MATLAB/mtex-6.0.0")
startup_mtex

%% Set up

% % crystal symmetry %% Turns out this is not necessary if you don't
% specify CS!! (Charlie and I realized this on 5/15/2025)

% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');


%This sets the format of figures to export
exportimagetype='image';
extension='.pdf';

% Set path & file names
pname = "/Users/cljd/Documents/SourceTree/Kamaehu_evo_2025/Other_Data/EBSD/AC1_ctf";
file = "KA69_c400 - EBSD Data.ctf";

filename = erase(file, " - EBSD Data.ctf");
fname = fullfile(pname, file);


% This turns off the X Y labels on pole figures
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations); 

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');

 % Charlie's EBSD correction
    % rotate everything 180 degrees around x to get the correct map
    ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
    % then rotate only the orientations 180 degrees around z
    ebsd = rotate(ebsd,rotation('axis',zvector,'angle',180*degree),'keepXY');

% Set up crystal shape for plotting
cS = crystalShape.olivine(ebsd('Forsterite').CS);

%% Check that map looks correct
fig=figure('Color','white');
[~,mP]=plot(ebsd);

set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k', 'Interpreter', 'latex');

savefig(fig,fullfile(pname, filename + "_RawPhaseMap.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_RawPhaseMap"+ extension), ...
    'ContentType', exportimagetype,'BackgroundColor','white');

%% Calculate grains
%Standard grain separation (10 degrees)
[grains, ebsd.grainId] = calcGrains(ebsd,'alpha',2.2,'angle',10*degree,'minPixel',10);

% Kuwahara 5 neighbor smoothing
ebsd(grains(grains.grainSize<2)).phase=0; %removes misindexed pixels 
F = KuwaharaFilter; 
F.numNeighbours = 5;
ebsd= smooth(ebsd('indexed'), F, 'fill',grains);

% Extra stringent grain definition (non standard)
subgrain_angle=0.3; %Subgrain definition angle
grain_div_angle=3; % Grain definition angle
len_filter = 2; % This is PLOT only subgrain boundaries with >2 length

[grains,ebsd.grainId]= calcGrains(ebsd,'alpha',2.2,'angle',[subgrain_angle*degree,grain_div_angle*degree],'minPixel',10); % Want this big, so it pulls out each grain separatly. 

% Smooth single phase grain boundaries
foGrains = grains('fo');

iter=[1,5,8];
for Z=1:length(iter)
    foGrains=smooth(foGrains, iter(Z));
end
foGrains=foGrains(foGrains.area >5000);

% filter inner boundaries to plot only long ones
gb = foGrains.innerBoundary;
long_gb = gb(gb.segLength > len_filter);

%% This is to save data and make the grain plot to match with Raman

numGrains = numel(foGrains);
% Set the variables we want to export
grainID            = foGrains.id;
GOS                = foGrains.GOS ./ degree;
grainSize          = foGrains.grainSize;
area               = foGrains.area;
subBoundaryLength  = foGrains.subBoundaryLength;
equivRadius        = foGrains.equivalentRadius;

% Create empty mats to store the GB lengths.
Tot_length  = zeros(numGrains,1);
prop_length = zeros(numGrains,1);

% Calculate the GB lengths per grain
for k = 1:numGrains
    gb_C_m = foGrains(k).innerBoundary;
    Tot_length(k) = sum(gb_C_m.segLength);
    prop_length(k) = Tot_length(k) / sqrt(grainSize(k));
end

placeholderCols = zeros(size(grainID)); % This is because right now it isn't calculating tilt and twist.
subgrain_angle_export  = subgrain_angle*(ones(size(grainID)));  
grain_div_angle_export = grain_div_angle*(ones(size(grainID)));  

Filename=repmat(string(filename), numel(grainID), 1);

% Create table
resultsTable = table(Filename,grainID, GOS, grainSize, area, subBoundaryLength, ...
    equivRadius, prop_length, Tot_length,...
    placeholderCols, placeholderCols, placeholderCols, ...
    placeholderCols, placeholderCols, ...
    subgrain_angle_export, grain_div_angle_export, ...
    'VariableNames', {'Filename','grainID',	'GOS' ,'Grain Size (pixels)', 'Grain Size (um2)', ...
    'subBoundaryLength (um)', 'equivalentRadius (um)','GB Length/Sqrt Size', 'Total GB length', ...
    'Tilt length',	'Twist length',	'Perc Tilt',	'Perc Twist', 'Perc unclassified', 'subgrain_angle', 'grain_div_angle'});


% Export to Excel
writetable(resultsTable, fullfile(pname,filename+'_grain_results.xlsx'));


%% Now we plot to match grains
fig=figure;
[~,mP]=plot(foGrains);
set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');

hold on
grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 12,'Color','black');
end
legend off
hold off
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k', 'Interpreter', 'latex');

savefig(fig,fullfile(pname, filename + "_Grains.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_Grains"+ extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');

%% All options past this are optional. Run sections individually.

% You can call the crystal symmetry either by using ebsd('fo').CS or by
% specifying a symmetry like so:
csEn = crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite  Opx AV77');
csFo = crystalSymmetry('mmm', [4.7560 10.2070 5.9800], 'mineral', 'Forsterite');

%% --- Basic maps ---

%% Nice band contrast map
fig=figure;
[~,mP]=plot(ebsd,ebsd.bc);
colormap gray
cb = mtexColorbar;
cb.Color = 'k';


set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');

hold on
plot(foGrains.boundary,'linewidth',1)
plot(long_gb, 'lineWidth', 0.5, 'lineColor', 'k');

grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 12);
end
legend off
hold off
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k', 'Interpreter', 'latex');

savefig(fig,fullfile(pname, filename + "_BandContrast.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_BandContrast"+ extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');

%% --- Internal distortion ---

%% Axis-angle colouring
% Good for visualising subtle distortion 
% Difficult to interpret the colours intuitively

% Set up colour key
AxisAngleKey = axisAngleColorKey(ebsd('fo'));
AxisAngleKey.oriRef = grains(ebsd('fo').grainId).meanOrientation;
AxisAngleKey.maxAngle = 3*degree;

color = AxisAngleKey.orientation2color(ebsd('fo').orientations);

fig=figure;
[~,mP]=plot(ebsd('fo'),color);
set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');


hold on
plot(foGrains.boundary,'linewidth',1)
plot(long_gb, 'lineWidth', 0.5, 'lineColor', 'k');

grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 12);
end
legend off
hold off
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k', 'Interpreter', 'latex');

savefig(fig,fullfile(pname, filename + "_Axisangle.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_Axisangle" + extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');

%% Angular intragranular misorientation in degrees (GROD)

% Calculate misorientation of each pixel relative to the average
% orientation of its grain
grod = ebsd('fo').calcGROD(grains);

fig=figure;
[~,mP]=plot(ebsd('fo'),grod.angle./degree);

set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');

mtexColorbar('title','Grain Reference Orientation Deviation - GROD angle (°)','Color','black')
setColorRange([0 3])
mtexColorMap(flipud(magma))%plasma%LaboTeX%parula
hold on
plot(foGrains.boundary,'linewidth',1)
plot(long_gb, 'lineWidth', 0.5, 'lineColor', 'k');

grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 12);
end

legend off
hold off
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k', 'Interpreter', 'latex');

savefig(fig,fullfile(pname, filename + "_GROD.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_GROD"+ extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');

%% IPF map
% Set up colour key
ipfKey = ipfHSVKey(ebsd('fo'));
ipfKey.inversePoleFigureDirection = orientation(grains(ebsd('fo').grainId).meanOrientation)*ipfKey.whiteCenter;
deg_rot=0;
ipfKey.colorPostRotation = rotation.byAxisAngle(zvector,deg_rot*degree);
ipfKey.maxAngle = 3*degree;
color = ipfKey.orientation2color(ebsd('fo').orientations);


fig=figure;
[~,mP]=plot(ebsd('fo'),color);

set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')
[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
    mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
    deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');


hold on
plot(foGrains.boundary,'linewidth',1)
plot(long_gb, 'lineWidth', 0.5, 'lineColor', 'k');
plot(grains('fo'),0.3*cS,'linewidth',0.1,'colored')

axis image off;

grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 12);
end
legend off
hold off

set(gca,'Color','white','Box','off') 
axis off
h = findobj(gca,'-property','DisplayName');
for k = 1:numel(h)
    h(k).DisplayName = erase(h(k).DisplayName,'$');

end

lgd = legend;
lgd.Interpreter = 'none';
lgd.Position = [0.7 0.3 0.05 0.1];
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
lgd.FontName = 'Times New Roman';
lgd.FontSize = 12;
ax1=gca;

nextAxis(2)
ipfKey_global = ipfHSVKey(ebsd('Forsterite').CS);
ipfKey_global.colorPostRotation=rotation.byAxisAngle(zvector,deg_rot*degree);
plot(ipfKey_global,'Title','none','ShowCoordinates','off');
axis image off;

ax2 = gca;

ax2.Units = 'normalized';
ax2.Position = [0.1 0.35 0.15 0.15];
ax2.FontName = 'Times New Roman';
ax2.FontSize = 12;
ax2.Color = 'none';
ax2.XColor = 'k';
ax2.YColor = 'k';
ax2.Title.Color='k';
uistack(ax1,'top')


savefig(fig,fullfile(pname, filename + "_IPF.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_IPF" + extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');

%% Basic crystal shape

fig=figure;
plot(cS,'colored')
ax=gca;
ax.Position=ax.Position.*[0.8 10.5 0.25 0.25];
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'black';
exportgraphics(fig,  fullfile(pname, "global_cS" + extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');
%% Basic KAM map (to visualise subgrain boundaries)

% Tell mtex that your pixels are in a grid
ebsd = ebsd.gridify;

kam = ebsd.KAM / degree;
fig=figure;
[~,mP]=plot(ebsd,kam);
mtexColorbar('title','Kernel Average Misorientation (°)','Color','k')
setColorRange([0,3]) % can fiddle with this or remove it
mtexColorMap LaboTeX
set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');

hold on
plot(foGrains.boundary,'linewidth',1)
plot(long_gb, 'lineWidth', 0.05, 'lineColor', 'gray');

grainlabels4plot=foGrains;
for j = 1:length(grainlabels4plot)
    % Calculate the center of the grain for labelling purposes
    centerX = mean(grainlabels4plot(j).boundary.x);
    centerY = mean(grainlabels4plot(j).boundary.y);

    % Get the grain ID and GOS value
    grainId = grainlabels4plot(j).id;
    GOS = grainlabels4plot(j).GOS ./ degree;

    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('%d', grainId), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 12);
end

legend off
hold off
set(gca,'Color','white','Box','off') 
axis off
lgd = legend;
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k', 'Interpreter', 'latex');
savefig(fig,fullfile(pname, filename + "_KAM.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_KAM"+ extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');



%% --- CRYSTALSHAPES ---

% To design a crystalShape for your mineral I recommend 
% midat.org to find likely habits
% https://www.smorf.nl/draw.php to design a shape (use 'crystallographic'
% setting)

% Some minerals (e.g. olivine) have shapes pre-loaded into mtex


% define the crystal shape of Forsterite and store it in the variable cS
cS = crystalShape.olivine(ebsd('Forsterite').CS);

%plot(cS,'colored')

fig=figure;
[~,mP]=plot(ebsd('fo'),'facecolor','lightgrey');
set(fig,'DefaultAxesFontName','Times New Roman')
set(fig,'DefaultTextFontName','Times New Roman')

[ mP.micronBar.visible, mP.micronBar.length, mP.micronBar.lineColor, ...
  mP.micronBar.scanUnit, mP.micronBar.backgroundColor ] = ...
  deal('on', 100, 'black', 'um', 'none'); 

set(mP.micronBar.txt,'FontSize',20,'FontName','Times New Roman','Color','black');

set(fig,'Renderer','painters');  % keep vector pipeline

hold on
plot(foGrains.boundary,'lineWidth',1)
plot(grains('fo'),0.4*cS,'linewidth',0.1,'colored')

grainlabels4plot=foGrains;
legend off
hold off

set(gca,'Color','white','Box','off') 
axis off
h = findobj(gca,'-property','DisplayName');
for k = 1:numel(h)
    h(k).DisplayName = erase(h(k).DisplayName,'$');

end

lgd = legend;
lgd.Interpreter = 'none';
lgd.Position = [0.85 0.1 0.05 0.1];
lgd.Color = 'none';
lgd.TextColor = 'black';
lgd.EdgeColor = 'none';
lgd.FontName = 'Times New Roman';
lgd.FontSize = 12;
savefig(fig,fullfile(pname, filename + "_Crystalshapes.fig"))
exportgraphics(fig,  fullfile(pname, filename + "_crystalshapes"+ extension), ...
    'ContentType', exportimagetype,'BackgroundColor','none');


% %% Euler map
% 
% colorKey = BungeColorKey(ebsd('fo').CS);
% 
% plot(ebsd('fo'),colorKey.orientation2color(ebsd('fo').orientations),'micronbar','off')
% fig=gcf;
% exportgraphics(fig,  fullfile(pname, filename + "_Eulercolor"+ extension), 'ContentType', exportimagetype);
% 
% %% IPF map
% 
% % make the colour key
% ipfKey_fo = ipfColorKey(ebsd('fo').CS);
% % specify reference direction (e.g., xvector, yvector, zvector)
% ipfKey_fo.inversePoleFigureDirection = zvector; 
% 
% 
% plot(ebsd('fo'),ipfKey_fo.orientation2color(ebsd('fo').orientations),'micronbar','off')
% fig=gcf;
% exportgraphics(fig,  fullfile(pname, filename + "_IPFcolor"+ extension), 'ContentType', exportimagetype);
% 





% %% --- Pole figures ---
% % For all pole figures: 'hkl' for poles to planes
% %                       'uvw' for crystallographic directions
% 
% %                       'lower', 'upper' or 'complete' for hemispheres
% 
% %                       'noSymmetry' to remove antipodal symmetry
% 
% % Can choose to specify a set of planes or axes like so:
% axes = Miller({1,0,0},{0,1,0},{0,0,1},'uvw',csFo);
% planes = Miller({1,0,0},{0,1,0},{2,1,0},{1,1,1},'hkl',csFo);
% 
% %% Plain pole figures (point per pixel)
% 
% %plotPDF(ebsd('fo').orientations,Miller({1,0,0},{0,1,0},{0,0,1},'hkl',csFo),'lower','MarkerColor','k','MarkerSize',2)
% 
% plotPDF(ebsd('fo').orientations,axes,'lower','MarkerColor','k','MarkerSize',2)
% 
% %% Plain pole figures (point per grain)
% plotPDF(grains('fo').meanOrientation,axes,'lower','MarkerColor','k','MarkerSize',2)
% 
% %% Axes coloured by ipfkey
% plotPDF(ebsd('fo').orientations,ipfKey_fo.orientation2color(ebsd('fo').orientations),Miller({1,0,0},{0,1,0},{0,0,1},'uvw',csFo),'lower','MarkerSize',2)
% 
% %% Planes coloured by ipfkey
% plotPDF(ebsd('fo').orientations,ipfKey_fo.orientation2color(ebsd('fo').orientations),Miller({1,0,0},{0,1,0},{0,0,1},'hkl',csFo),'lower','MarkerSize',2)
% 
% %% Smoothed pole figure
% 
% % Smoothed
% plotPDF(ebsd('fo').orientations,Miller({1,0,0},{0,1,0},{0,0,1},'hkl',csFo),'lower','smooth','halfwidth',5*degree)
% mtexColorMap LaboTeX
% setColorRange('equal')
% mtexColorbar('location','eastoutside','title','multiples of uniform distribution')