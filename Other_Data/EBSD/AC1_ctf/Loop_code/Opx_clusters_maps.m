% Opx cluster loop: maps
 
% Charlie Gordon 20th Feb 2025
% Alex's data, collected on the Berkeley SEM in Jan 2025

% Desgined to run in mtex 6.0.0 (uses updated grain calculation code)

close all
clear all

% Start up mtex
addpath("C:\Users\charl\OneDrive\Documents\MATLAB\mtex-6.0.0\mtex-6.0.0")
startup_mtex

%% Set up
% csEn = crystalSymmetry('mmm', [18 8.8 5.2], 'mineral', 'Enstatite  Opx AV77');

% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Set path & file names

%pname = 'C:\Users\charl\Dropbox\SEM\EBSD data\0_not_my_data\20250128_Alex_denoised_ctfs';

pname = "P:\WORK-GENERAL\POSTDOC-UCB\BERKELEY-VIBE\Documents\Projects\" + ...
    "PrePostCaldera_Kil2025\EBSD\Kilauea_EBSD_05222025\ctf";

fname = dir(fullfile(pname, '*.ctf'));



%% Run loop

% Set up loop
for k = 1:length(fname)
    baseFileName = fname(k).name;
    fullFileName = fullfile(pname, baseFileName);

    % Extract the sample name ('name') from the extension
    [pathstr,name,ext] = fileparts(baseFileName);
   
    %Import the Data
    ebsd = EBSD.load(fullFileName,'interface','ctf',...
    'convertEuler2SpatialReferenceFrame');

    % Charlie's EBSD correction
    % rotate everything 180 degrees around x to get the correct map
    ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
    % then rotate only the orientations 180 degrees around z
    ebsd = rotate(ebsd,rotation('axis',zvector,'angle',180*degree),'keepXY');

    % Simple Enstatite crystalShape
    Fo_cS = crystalShape.olivine(ebsd('Forsterite').CS);

    % Calculate grains
    [grains, ebsd.grainId] = calcGrains(ebsd,'alpha',2.2,'angle',10*degree,'minPixel',5);
    [grains, ebsd.grainId] = calcGrains(ebsd,'alpha',2.2,'angle',3*degree,'minPixel',5);


    % --- Approach using grain boundary segments ---
    % This calculates the grain to grain misorientation (for grains touching each other)
    % % Extract grain boundaries of interest
    % gB = grains.boundary;
    % gB_Fo_Fo = grains.boundary('Forsterite','Forsterite');
    % 
    % % Calculate the grain boundary misorientations
    % diso = gB_Fo_Fo.misorientation;
%
    % % COMBINED FIGURE
    % close all
    % 
    % % Layout is rows then columns
    % newMtexFigure('figSize','large','layout',[1,2])
    % 
    % 
    % % gB coloured by misorientation axis (with crystalShapes)
    % 
    % mis_axes = gB_Fo_Fo.misorientation.axis;
    % 
    % colorKey = HSVDirectionKey(csEn);
    % color = colorKey.direction2color(mis_axes);
    % 
    % plot(ebsd('fo'),'facecolor','lightgrey','micronbar','off')
    % hold on
    % plot(gB_Fo_Fo,color,'linewidth',4)
    % plot(grains('fo'),0.5*Fo_cS,'linewidth',0.1,'colored')
    % 
    % legend off
    % hold off
    % title('gB coloured by mis axes','FontSize',10)
    % 
    % % Grain boundary misorientation axes on the colour key
    % 
    % nextAxis
    % plot(colorKey)
    % hold on
    % plot(mis_axes,'MarkerSize',1.5,'MarkerColor','k')
    % hold off
    % title('gB misorientation axes','FontSize',10)
    % 
    % 


 % Plot of grod with other phases in grey
    
    figure

    grod = ebsd('fo').calcGROD(grains);

    plot(ebsd('indexed'),'facecolor','darkgrey','micronbar','off')
    hold on
    plot(ebsd('fo'),grod.angle./degree,'micronbar','off')
    mtexColorbar('title',{'GROD angle in degrees (ol only)'})
    setColorRange([0 5])
    mtexColorMap LaboTeX

    plot(grains.boundary,'lineWidth',1)
    legend off
    hold off
    title('GROD angle','FontSize',10)

     
    % --- Save ---
    % Save figure
    saveas(gcf,[name, '_grod.png'])
 
    
    close all


 end


%%

 % Grain boundaries coloured by angle
    
    % nextAxis
    % plot(ebsd('fo'),'facecolor','lightgrey','micronbar','off')
    % hold on
    % plot(gB_Fo_Fo,gB_Fo_Fo.misorientation.angle./degree,'linewidth',4)
    % mtexColorbar('title','misorientation angle (°)')
    % legend off
    % hold off
    % title('gB coloured by mis angle','FontSize',10)
    % 

%% 

%plot(ebsd)