%% Calculate the cLh, cLa


%% Inputs

%%%Read the rest of the data from json
filename = 'user_input.json';
fid = fopen(filename, 'r');
rawData = fread(fid, '*char')';  % Read the file as a string
fclose(fid);

% Decode the JSON string into a MATLAB structure
data = jsondecode(rawData);

% Read all input variables from the string
Mach =data.Mach;
RedFreq=data.RedFreq;  %%defined as wb/U
b=data.ref_semichord;                      %Gust modeshape input is a column vector with dimensions: # of aero boxes on gust vanes only x 1

%%Add paths
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'DLM_helperfunctions'));
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'Spline_helperfunctions'));

%%%Get some data from the mesher:
[ProjectedBoxArea,BoxArea,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection] = DLM_Mesh_Updated();

%% Define modes
[AIC] = DLM_Find_AICs(Mach, RedFreq/b,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection); 

verticalmode_at_downwash=ones(264,1); %upward motion by 1m.
verticalmode_at_forcepoints=ones(264,1);

rotationmode_at_downwash=BoxPointDownwash(:,1); % deflection due to rotation by 1 rad around x=0 ( 1 rad pitch)
rotationmode_at_downwash=BoxSendPoint3(:,1);

zout_div_dx_dw_verticalmode=dTdx_Stickmodel(BoxPointDownwash(:,1),BoxPointDownwash(:,2),verticalmode_at_downwash,BoxPointDownwash(:,1),BoxPointDownwash(:,2));
zout_div_dx_dw_rotationmode=dTdx_Stickmodel(BoxPointDownwash(:,1),BoxPointDownwash(:,2),rotationmode_at_downwash,BoxPointDownwash(:,1),BoxPointDownwash(:,2));

%% Work done on mode 1 (ie: the lift direction) by mode 2 (one radian of angle of pitch):

Q12=verticalmode_at_forcepoints'*diag(ProjectedBoxArea)*inv(AIC)*(1i*RedFreq/b*rotationmode_at_downwash+zout_div_dx_dw_rotationmode); 

Cl_alpha=Q12/sum(ProjectedBoxArea);

%% Work done on mode 1 (ie: the lift direction) by mode 1

Q11=verticalmode_at_forcepoints'*diag(ProjectedBoxArea)*inv(AIC)*(1i*RedFreq/b*verticalmode_at_downwash+zout_div_dx_dw_verticalmode); 

Cl_h=Q11/sum(ProjectedBoxArea);

disp(Cl_alpha)
disp(Cl_h)