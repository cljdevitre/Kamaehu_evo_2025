clear all 
close all
addpath('C:\Users\penny\Box\Berkeley_new\EBSD\EBSD_Codes\mtex-5.11.1');
startup_mtex;
subgrain_angle=0.3;
grain_div_angle=4;


% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

AAoutput={}


%% First read from the spreadsheet

% Load the spreadsheet
data = readtable('Merged_GOS_parameters_edited.xlsx');

% Filter rows where 'Slip system' is 'Y'
filtered_data = data(strcmp(data.Slip_system, 'Y'), :);

% Loop through each remaining row
for i = 51:height(filtered_data)
    disp('Looping')
    disp(i)

    % Set fname to the value in the column 'Filename'
    f_to_save=filtered_data.Filename{i}
    fname = filtered_data.Filename{i};
    fname = strrep(fname, '/', '');
    fname = strrep(fname, '\', '');
    fname = strrep(fname, 'EBSD', ''); % Charlotte and Christine may find this causes issues
    fname = regexprep(fname, 'Map.*', '')
    % Use regexprep to remove all spaces
    fname = regexprep(fname, '\s+', '');
    day = filtered_data.Day(i);
    switch day
    case 1
        pname = 'C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\OldEBSDProjects_beforeReexport\Data_ML_EBSD\Day1';
    case 2
        pname = 'C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\OldEBSDProjects_beforeReexport\Day2_MLB';
    case 3
        pname =  'C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\OldEBSDProjects_beforeReexport\Day three EBSD\EBSD ML Day 3\ctf_files_day3';
    case 4
        pname =  'C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\day four ebsd\Day 4 ebsd\ctf';
    case 5
        pname='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD\Day five Christina + ML\Ori1';

    otherwise
        error('Invalid day value. Please specify day as 1-5.');
    end

    
% Get a list of all files in the directory
files = dir(pname);


% Loop through each file and check if the name contains the string 'fname'
file = '';
for j = 1:length(files)
    if contains(files(j).name, fname)
        file = fullfile(pname, files(j).name);
        break;
    end
end








%%
    if contains(file, 'ctf')
        ebsd = EBSD.load(file, 'convertEuler2SpatialReferenceFrame', 'interface', 'ctf');
        % rotate everything 180 degrees around x to get the correct map
        ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
        % then rotate only the orientations 180 degrees around z
        ebsd = rotate(ebsd,rotation('axis',zvector,'angle',180*degree),'keepXY');
    
    else
        
        ebsd = EBSD.load(fname,'interface','h5oina', 'convertEuler2SpatialReferenceFrame');
% rotate everything 180 degrees around x to get the correct map
    ebsd = rotate(ebsd,rotation('axis',xvector,'angle',180*degree));
    ebsd = rotate(ebsd,rotation('axis',yvector,'angle',180*degree),'keepXY');

    end
    ebsd.unitCell = ebsd.unitCell;

%% Enter strike of boundaries with a given WBV
% in our coordinate system, x is to the right, y is up, z is out of the
% page
% 1, 0, 0 would be east. 
% 0, 1, 0 would be north. 
% 100 is actually blue in Aztec crystal now. 
    % Set 'greenstrike' based on the value of 'Blue direction'
   blueDirection = filtered_data.Blue_direction{i};
    switch blueDirection
        case 'N'
            greenstrike = vector3d(0, 1, 0);
        case 'NE'
            greenstrike = vector3d(1, 1, 0);
        case 'E'
            greenstrike = vector3d(1, 0, 0);
        case 'SE'
            greenstrike = vector3d(1, -1, 0);
        otherwise
            error('Invalid Blue direction value.');
    end 
    % Read the 'Red direction' column
    redDirection = filtered_data.Red_direction{i};
    
    % Set 'redstrike' based on the value of 'Red direction'
    if isempty(redDirection)
        redstrike = NaN;
    else
        switch redDirection
            case 'N'
                redstrike = vector3d(0, 1, 0);
            case 'NE'
                redstrike = vector3d(1, 1, 0);
            case 'E'
                redstrike = vector3d(1, 0, 0);
            case 'SE'
                redstrike = vector3d(1, -1, 0);
            otherwise
                error('Invalid Red direction value.');
        end
    end

%% Data Filtering
[grains, ebsd.grainId, ebsd.mis2mean]= calcGrains(ebsd,'angle',10*degree);
ebsd(grains(grains.grainSize<2)).phase=0; %removes misindexed pixels 
% Apply 5 neighbour Kuwahara filter
F = KuwaharaFilter; 
F.numNeighbours = 5;
ebsd= smooth(ebsd('indexed'), F);
%% Subgrain reconstruction
[grains,ebsd.grainId,ebsd.mis2mean]= calcGrains(ebsd,'angle',[subgrain_angle*degree,grain_div_angle*degree]);


% loops over smoothing multiple times.
iter=[1,5,8];
for Z=1:length(iter)
    grains=smooth(grains, iter(Z));
end
grains=grains(grains.area >1000); % grains(grains.grainSize > max(grains.grainSize/30))  % grains(grains.area > 100); % Only use grains >100 um^2 in size (so 10X 10 um) % 
grains=grains('Fo'); %Construct grains variable that is restricted to only the phase forsterite

%%
GOS_target = filtered_data.GOS(i);
    tolerance = 0.009; % Define an appropriate tolerance value
    min_difference = Inf; % Initialize min_difference to a large value
    closest_grain = -1; % Initialize closest_grain to an invalid value

    for k = 1:length(grains) % Inner loop
        disp('looping through grains to find right GOS');
        GOS = grains(k).GOS / degree;
        display(GOS);
        display(GOS_target);
        GOS_diff = abs(GOS - GOS_target);
        
        if GOS_diff < min_difference
            min_difference = GOS_diff;
            closest_grain = k;
            display(GOS_diff);
            if GOS_diff < tolerance
                break; % Exit the inner loop if the difference is within the tolerance
            end
        end
    end

    % Check if we found a closest grain
    if closest_grain == -1
        disp('No grain found within the specified tolerance.');
        continue; % Move to the next i iteration
    else
        disp(['Closest grain found at index: ', num2str(closest_grain)]);
        disp(['Target GOS: ', num2str(GOS_target)]);
        disp(['GOS of target grain: ', num2str(GOS)]);
        
        % Check if the closest GOS is still outside tolerance
        if min_difference >= tolerance
            disp(['Closest GOS difference is too large: ', num2str(min_difference)]);
            continue; % Move to the next i iteration
        end
    end
    


%% Visually inspect EBSD data, grains and grain boudnaries
figure(1)

hold on

    currentGrainIdsz = grains(k).id;  % IDs of pixels in the current grain
    

     % Combine grain ID with filename

    gb_C = grains('id',currentGrainIdsz).innerBoundary;
    cs = ebsd('Fo').CS
    mAng = gb_C.misorientation.angle;

    % Misorientation axis in crystal coordinates
    gB_axesC = gb_C.misorientation.axis;
    %axis in specimen coords
    o = ebsd('id',gb_C.ebsdId).orientations;
    ROT = axis(o(:,1),o(:,2))
    ROT.antipodal=1;
    %Calculate grain boundary directions in specimin coordinates
    GBD=gb_C.direction;

    %
     oM = ipfHSVKey(ebsd(ebsd.grainId == currentGrainIdsz));
    oM.inversePoleFigureDirection = orientation(grains('id',currentGrainIdsz).meanOrientation) * oM.whiteCenter;
    oM.colorPostRotation = reflection(yvector);
    oM.maxAngle = 2*degree;
    color=oM.orientation2color(ebsd(ebsd.grainId == currentGrainIdsz).orientations);
    hold on
    [~,mP]=plot(ebsd(ebsd.grainId == currentGrainIdsz), color);
   

   

% plot
% plot(grains('id',currentGrainIdsz).boundary)
% hold on
% ckc = HSVDirectionKey(gB_axesC.CS)
% plot(gb_C,ckc.direction2color(axC(cond)))
% hold off


%Calculate the direction of the misorientation axis in specimin coordinates. 
% I think this is the issue - got to here on friday night. Doesnt seem to https://github.com/mtex-toolbox/mtex/discussions/2179


%%
figure(33)
plot(ebsd(ebsd.grainId == currentGrainIdsz), ebsd(ebsd.grainId == currentGrainIdsz).bc);
hold on
plot(gb_C, 'linecolor','k','linewidth',1.5) % This is all boundaries.
%% Remove high angle grain boundaries, and very short boundaries
condnoise=gb_C.misorientation.angle./degree>subgrain_angle;
condbig=gb_C.misorientation.angle./degree<grain_div_angle;
% Only include boundaries longer than 2 scan units. 
condlength=gb_C.segLength>0; 

%% Considering boundary dip - shallowly dipping boundaries are unlikely to be observed
% Calculate vector perpendicuular to misorientation axis and grain boundary direction in specimin coordinates. This is the normal to a tilt boundary
tiltnormal=cross(ROT, GBD);
% Define vertical vector in specimin coordinates
z=vector3d(0, 0, 1);
% Condition tilt 1: Only consider tilt boundaries that dip steeper than 15
% degree wrt. the vertical. 90 degrees would be the perfect tilt boundary
condtilt1=angle(tiltnormal,z, 'antipodal')>15*degree;
% Condition twist 1 : only consider twist boundaries that dip steeper than
% 15 degree wrt. the vertical. 0 would be looking down the twist axis.
condtwist1=angle(ROT,z, 'antipodal')>15*degree;
% Condition twist 2 : Ideal twist boundary has boundary strike perpendicular to rotation axis. Here we allow 15 degrees leniancy. 
condtwist2=angle(ROT, GBD, 'antipodal')>75*degree;
%% Classify boundaries as tilt or twist
% If a boundary meets both twist criteria, classify as a twist boundary. Include angular, noise and length critera. 
TwistLog=condtwist1&condtwist2;
TwistBoundaries=TwistLog&condnoise&condbig&condlength;
% Classify all boundaries that are not twist as tilt boundaries, Include angular, noise and length critera. 
TiltLog=~TwistLog;
TiltBoundaries=TiltLog&condnoise&condbig&condlength&condtilt1;
%% Plot figure showing boundaries classifed as tilt and twist
figure(3)
mtexFig = newMtexFigure;


[~,mP]=plot(ebsd(ebsd.grainId == currentGrainIdsz), color);



hold on
plot(gb_C(TiltBoundaries), 'linecolor','r','linewidth',2,'DisplayName', 'Tilt');
plot(gb_C(TwistBoundaries), 'linecolor','y','linewidth',2,'DisplayName', 'Twist');
mP.micronBar.visible = 'off';
legendHandle = legend('show');
set(legendHandle, 'FontSize', 30);
%% MAke a folder for slip system
pname_pic='C:\Users\penny\Box\Berkeley_new\Raman\MaunaLoa\EBSD'
    folder_tt = fullfile(pname_pic, 'SlipSystemPics');
    if ~exist(folder_tt, 'dir')
        mkdir(folder_mis);
    end
savePath = fullfile(folder_tt, [strrep(fname, '.ctf', ''), '_TiltVsTwist.png']);
saveas(gcf, savePath)

%% Classifying tilt boundaries with [100] WBV directions by fabric type, and progressively plot onto this figure
figure(4)
mtexFig = newMtexFigure;
 [~,mP]=plot(ebsd(ebsd.grainId == currentGrainIdsz), color);


mP.micronBar.visible = 'off'
hold on
plot(gb_C)
%

% Isolating bounadries with different WBV directions (leniancy of 40 degrees to the entered direction)
condgreenstrike=angle(greenstrike, GBD, 'antipodal')<40*degree;

% Defining tilt boundaries with A type fabrics 
%A type fabrics have misorientation axis in crystal coordinates within 20 degrees of the [001] axis
cond2=angle(gB_axesC, Miller(0,0,1,cs))<20*degree;
% A tilt boundaries meet this crystallographic criteria and have [100] burgers vectors
potentialATiltBoundaries = gb_C(cond2&condgreenstrike&TiltBoundaries)
%If any A tilt boundaries are found, plot on figure
if length(potentialATiltBoundaries)>1
    % IT is cond2 that breaks this  -But cond2 works on its own.  
    combinedCondition1 = cond2 & condgreenstrike;
combinedCondition2 = combinedCondition1 & TiltBoundaries;
plot(potentialATiltBoundaries, 'linecolor','cyan','linewidth',4, 'DisplayName', 'Atype');
end



%% Defining tilt boundaries with E type fabrics 
%E type fabrics have misorientation axis in crystal coordinates within 20 degrees of the [010] axis
cond3=angle(gB_axesC, Miller(0,1,0,cs))<20*degree;
% E tilt boundaries meet this crystallographic criteria and have [100] burgers vectors
potentialETiltBoundaries = gb_C(cond3&condgreenstrike&TiltBoundaries);
%If any E tilt boundaries are found, plot on figure
if length(potentialETiltBoundaries)>1;
plot(potentialETiltBoundaries, 'linecolor','blue','linewidth',4 , 'DisplayName', 'Etype');;
end
%% Defining tilt boundaries with D type fabrics 
%D type fabrics have misorientation axis in crystal coordinates >70 degrees from [100], and more than 20 degrees away 
%from [010] and [001] (else would be A or E. 
cond4=angle(gb_C.misorientation.axis, Miller(1,0,0,cs))>70*degree;
cond5=angle(gb_C.misorientation.axis, Miller(0,1,0,cs))>20*degree;
cond6=angle(gb_C.misorientation.axis, Miller(0,0,1,cs))>20*degree;
% D tilt boundaries meet this crystallographic criteria and have [100] burgers vectors
potentialDTiltBoundaries = gb_C(cond4&cond5&cond6&condgreenstrike&TiltBoundaries);
%If any D tilt boundaries are found, plot on figure
if length(potentialDTiltBoundaries)>1;
plot(potentialDTiltBoundaries, 'linecolor','yellow','linewidth',4, 'DisplayName', 'Dtype');
end
legendHandle = legend('show');
set(legendHandle, 'FontSize', 30);
savePath = fullfile(folder_tt, [strrep(fname, '.ctf', ''), '_100WBVD_Boundaries.png']);
saveas(gcf, savePath)

%% Classifying tilt boundaries with [001] WBV directions by fabric type, and progressively plot onto this figure
if ~isnan(redstrike)
    condredstrike=angle(redstrike, GBD, 'antipodal')<40*degree;
    figure(5);
    mtexFig = newMtexFigure;
for z = 1:length(grains)

    currentGrainIdsz = grains(z).id;  % IDs of pixels in the current grain

     % Combine grain ID with filename

    currentGrainEBSDz = ebsd(ebsd.grainId == currentGrainIdsz);  % Subset EBSD data for current grain
    oM = ipfHSVKey(currentGrainEBSDz);
    oM.inversePoleFigureDirection = orientation(grains(z).meanOrientation) * oM.whiteCenter;
    oM.colorPostRotation = reflection(yvector);
    oM.maxAngle = 2*degree;
    color=oM.orientation2color(currentGrainEBSDz.orientations);
    hold on
    [~,mP]=plot(ebsd(ebsd.grainId == currentGrainIdsz), color);

end
hold on
plot(gb_C)
mP.micronBar.visible = 'off'
%% Defining tilt boundaries with B type fabrics 
%B type fabrics have misorientation axis in crystal coordinates within 20 degrees of the [100] axis
condB=angle(gb_C.misorientation.axis, Miller(1,0,0,cs))<20*degree;
% B tilt boundaries meet this crystallographic criteria and have [001] burgers vectors
potentialBTiltBoundaries = gb_C(condB&condredstrike&TiltBoundaries)
%If any B tilt boundaries are found, plot on figure
if length(potentialBTiltBoundaries)>1
plot(potentialBTiltBoundaries, 'linecolor','green','linewidth',4, 'DisplayName', 'Btype');
end
%% Defining tilt boundaries with C type fabrics 
%C type fabrics have misorientation axis in crystal coordinates within 20 degrees of the [010] axis
condC=angle(gb_C.misorientation.axis, Miller(0,1,0,cs))<20*degree;
% C tilt boundaries meet this crystallographic criteria and have [100] burgers vectors
potentialCTiltBoundaries = gb_C(condC&condredstrike&TiltBoundaries)
%If any C tilt boundaries are found, plot on figure
if length(potentialCTiltBoundaries)>1
plot(potentialCTiltBoundaries, 'linecolor','red','linewidth',4, 'DisplayName', 'Ctype');
end
legendHandle = legend('show');
set(legendHandle, 'FontSize', 30);
savePath = fullfile(folder_tt, [strrep(fname, '.ctf', ''), '_001WBVD_Boundaries.png']);
saveas(gcf, savePath)

end
%% Calculating Lengths of each boundary type
TiltLength=sum(gb_C(TiltBoundaries).segLength);
TwistLength=sum(gb_C(TwistBoundaries).segLength);
if length(potentialATiltBoundaries)>0;
Alength=sum(potentialATiltBoundaries.segLength);
else
    Alength=0;
end

if length(potentialDTiltBoundaries)>0;
Dlength=sum(potentialDTiltBoundaries.segLength);
else
    Dlength=0;
end

if length(potentialETiltBoundaries)>0;
Elength=sum(potentialETiltBoundaries.segLength);
else
    Elength=0;
end

if ~isnan(redstrike);
if length(potentialBTiltBoundaries)>0;
Blength=sum(potentialBTiltBoundaries.segLength);
else
    Blength=0;
end

if length(potentialCTiltBoundaries)>0    ;
Clength=sum(potentialCTiltBoundaries.segLength);
else
    Clength=0;
end

else
    Blength=0;
    Clength=0;
end
%% 

AAoutput=[AAoutput; {i, f_to_save, file, day, filtered_data.GOS(i), Alength Dlength Elength Blength Clength TiltLength TwistLength}]; 
   % if mod(i, 10) == 0
   %      Convert AAoutput to a table
   %      T = cell2table(AAoutput, 'VariableNames', {'File', 'Day', 'GOS', 'Alength', 'Dlength', 'Elength', 'Blength', 'Clength', 'TiltLength', 'TwistLength'});
   % 
   %      Write to Excel
   %      writetable(T, excelFile, 'Sheet', 'Sheet1', 'WriteMode', 'append', 'WriteVariableNames', false);
   % end
        %Store in excel to aggregate various boundaries
close all

varsToKeep = {'files', 'folder', 'pname', 'i', 'subgrain_angle', 'grain_div_angle', 'color_angle', 'AAoutput', 'filtered_data'};
allVars = who;% Get a list of all variables in the workspace


varsToClear = setdiff(allVars, varsToKeep); % Determine which variables to clear

clear(varsToClear{:}); % Clear the selected variables


end


