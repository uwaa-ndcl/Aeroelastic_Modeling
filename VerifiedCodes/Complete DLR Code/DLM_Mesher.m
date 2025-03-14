
% Lifting Surface by the Doublet Lattice Method
%
% EL November 7 2024
%
% Meshing: Creating the array of boxes based on panel geometry information.
%
%
% At this stage, geometry data is entered here directly into the source code.
% To make the code more professional, data should be read from a data file
% and users will not touch the code itself.
%
% Incoming flow - along the x axis
%
clear all
close all
clc
tic
im=1i; % the pure imaginary number
%
% Defining the geometry of trapezoidal panels
%
% Axes:
% x from nose to tail in the direction of the flow
% y to the right
% z up 
%
%
%   |
%   |
%   |
%   | U
%   |
%   |---------------------------------- y
%   |
%   |  *FL
%   |   |
%   |   |               *FR
%   |   |               |
%   |   |               *AR
%   |  *AL
%   |
%   x
%
% Each panel is defined by the coordinates of 4 corner points
% 1-2-3-4 
% from the forward left (FL) 
% to the forward right (FR)
% to the aft-right (AR) 
% to the aft-left (AL) 
% in the clockwise direction.
%
% The array data_panel stores the x,yz coordinates of points FL,FR,AR,AL
% for each panel plus the number of boxes in the chordwise and spanwise
% directions that will cover the panel.
% Each row corresponds to a different panel
% Each row includes: 
% xFL yFL zFL xFR yFR zFR xAR yAR zAR xAL yAL zAL nchord nspan
%
% The array PanelPoints stores the coordinates of the 4 points per
% Panel as follows:
% each column stores x,y,z coordinates of a point
% column1 - point 1, an so on to column 4 - point 4
% a third index identifies the Panel.
% And so PanelPoints(i,j,k) contains the i-th coordinate (x,y, or z)
% of the j-th point (FL, FR, AR, or AL) of the k-th Panel.
%
% ------------------------------------------------------------------
% In the following: a few test cases.
% All are commented out except for just one selected to run.
% ------------------------------------------------------------------
%
% Cornell Wing
%
% ' Cornell Wing' 
%
%data_panel= ...
%[-.22874 0.      0. .317674 .54648 0. .77521  .54648 0. .22874  0.      0. 6 6;
%  .45747 0.      0. 1.00395 .54648 0. 1.46142 .54648 0. .91494  0.      0. 6 6 ;
%  .31774 -.54648 0. -.22874 0.     0. .22874   0.    0. .77521  -.54648 0. 6 16 ;
% 1.00395 -.54648 0. .45747  0.     0. .91494   0.    0. 1.46142 -.54648 0. 6 16] ;
%refchord=1.0;
% Mach and reduced frequency
%data_aero=[ 0.0 2.6231 ];
%Mach=data_aero(1)
%RedFreq=data_aero(2)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------
% TEST CASE AR=2 Rectangular Wing
% EXAMPLE IN WL-TR-95-3022
%  FL(x,y,z)  FR(x,y,z)  AR(x,y,z)  AL(x,y,z) DivCh DivSp
%
% ' The AR=2 rectangular wing '
%
% data_panel=[0 -2 0 0 0 0 2 0 0 2 -2 0 10 10;
%             0  0 0 0 2 0 2 2 0 2  0 0 10 10];
% refchord=2.;

% Mach ReducedFrequency
% data_aero=[0.5 0.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------
% ------------------------------------------------------------------
% ' Wing+Flap Case '
% 
%           FL(x,y,z)  FR(x,y,z)  AR(x,y,z)  AL(x,y,z) DivCh DivSp
%
% data_panel=[0 -2 0 0 2 0 1 2 0 1 -2 0 15 15;
%            1 -2 0 1 2 0 2 2 0 2 -2 0 15 15];
% refchord=???
% Mach ReducedFrequency
% data_aero=[0.5 0.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------
%
%' AGARD WING-HORIZONTAL TAIL CONFIGURATION ' %          <-----------------Input
%      FL(x,y,z)  FR(x,y,z)  AR(x,y,z)  AL(x,y,z) DivCh DivSp
 % data_panel=[2.75 -1.  0. 0.    0. 0. 2.25  0. 0.   3.7 -1   0  9 10;
 %             0.0  0.0  0. 2.75  1. 0. 3.7   1. 0.   2.25 0.  0. 9 10;
 %             3.9  -1.  0. 2.7   0. 0. 4.    0. 0.   4.25 -1. 0. 9 10;
 %             2.7  0.0  0. 3.9   1. 0. 4.25  1. 0.   4.   0.  0. 9 10];
%  refchord=1.0 ;
%            Mach ReducedFrequency
% data_aero=[0.8 1.5];                     % <-------------------------------Input
% Mach=data_aero(1);    % Note: this is not needed for setting up the geometry
% RedFreq=data_aero(2); % Note: More than one reduced frequency should be allowed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ' SWEPT AR=2 WING IN WL-TR-95-3022 '
% % %           FL(x,y,z)  FR(x,y,z)  AR(x,y,z)  AL(x,y,z) DivCh DivSp
% data_panel=[0.81 -1 0 0 0 0 1.614 0 0 1.954 -1 0 10 10;
%             0 0 0 0.81 1 0 1.954 1 0 1.614 0 0 10 10];
% refchord=????;
% Mach ReducedFrequency
% data_aero=[0.8 1.4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ' TEST CASE 5%%%EXAMPLE IN AIAA 81-4157 '
%           FL(x,y,z)  FR(x,y,z)  AR(x,y,z)  AL(x,y,z) DivCh DivSp
% data_panel=[.318 -.546 0 -.229 0 0 .229 0 0 .775 -.546 0 8 12;
%     -.229 0 0 .318 .546 0 .775 .546 0 .229 0 0 8 12;
%     1.004 -.546 0 .457 0 0 .915 0 0 1.461 -.546 0 8 12;
%     .457 0 0 1.004 .546 0 1.461 .546 0 .915 0 0 8 12];
% refchord=?????
% Mach ReducedFrequency
% data_aero=[0.7 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------
%
% Forshing wing from Rowe's Paper
%
% ' 25 deg swept wing / flap ' 
%
%cos25=0.90630779;
% data_panel= ...
% [0.0 0.0 0.0 0.19119 0.41 0.0 0.61119 0.41 0.0 0.42 0.0 0.0 6 6 ;
%  0.19119 -0.41 0.0 0.0 0.0 0.0 0.42 0.0 0.0 0.61119 -0.41 0.0 6 6 ;
%  0.42 0.0 0.0 0.611190 0.41 0.0 0.79119 0.41 0.0 0.6 0.0 0.0 6 6 ;
%  0.61119 -0.41 0.0 0.42 0.0 0.0 0.6 0.0 0.0 0.791119 -0.41 0.0 6 6 ;
%  0.19119 0.41 0.0 0.19585 0.42 0.0 0.61585 0.42 0.0 0.61119 0.41 0.0 6 1 ;
%  .19585 -.42 0. 0.19119 -0.41 0. 0.61119 -0.41 0.0 0.61585 -0.42 0.0 6 1 ;
%  0.61119 0.41 0.0 0.61585 0.42 0.0 0.79585 0.42 0.0 0.79119 0.41 0.0 6 1;
%  .61585 -.42 .0 .61119 -0.41 .0 .79119 -.41 .0 .79585 -.42 .0 6 1 ;
%  .19585 .42 0. .41035 .88 0. .83035 .88 .0 .61585 0.42 .0  6 6 ;
%  .41035 -0.88 0.0 0.19585 -0.42 0.0 0.61585 -0.42 0.0 0.83035 -0.88 0.0 6 6;
%  0.61585 0.42 0.0 0.83035 0.88 0.0 1.01035 0.88 0.0 0.79585 0.42 0.0 6 6;
%  0.83035 -0.88 0. 0.61585 -0.42 .0 0.79585 -0.42 0.0 1.01035 -0.88 0. 6 6] ;
%refchord=1.0;
%data_aero=[ 0.0 2.6231 ];
%Mach=data_aero(1)
%RedFreq=data_aero(2)
% 
%
% ----------------------------------------------------------
% End of collection of test data cases 
% ----------------------------------------------------------



% ----------------------------------------------------------
% Beginning of reading data from .json file
% Comment this section out and enable one of the above test cases
% If you want to define the input directly in this Matlab code
% In addition, you will need to define NSmax, NCmax and plot_control
% variables in this code if you will not use the .json input file.
% ----------------------------------------------------------

% Load and decode the JSON file into a MATLAB structure
filename = 'user_input.json';
fid = fopen(filename, 'r');
rawData = fread(fid, '*char')';  % Read the file as a string
fclose(fid);

% Decode the JSON string into a MATLAB structure
data = jsondecode(rawData);

% Read all input variables from the string
data_panel= data.data_panel;
refcord =data.refchord;
Mach =data.Mach;
RedFreq=data.RedFreq;
NSmax=data.NSmax;
NCmax=data.NCmax;
plot_control=data.plot_control;
% ----------------------------------------------------------
% End of reading data from .json file 
% ----------------------------------------------------------


% Control of plotting (zeros or 1s):
% element 1 - plot panels & Aero Boxes
% element 2 - plot vortex edge points on boxes
% Element 3 - plot force points and downwash on boxes
% ---------------------
%plot_control=[1 0 1 ];   % This was previously inputted here. Now, it is
%read from the input .json file
% ---------------------
%
%
plotbox=plot_control(1);
plotVortex=plot_control(2);
plotforceanddownwash=plot_control(3);
%
% Angle of Attack
%
FlowDirection=[1. 0. 0.]';
%
% Geometry inputs:
%
nPanels=size(data_panel,1);% No. of Panels=length of column 1 of data_panel
%
% As already presented above:
% Each panel is defined by the coordinates of 4 corner points
% 1-2-3-4 from the forward left (FL) to the forward right (FR)
% to the aft-right (AR) to the aft-left (AL) in the clockwise
% direction.
%
% The array PanelPoints stores the coordinates of the 4 points per
% Panel as follows:
% each column stores x,y,z coordinates of a point
% column1 - point 1, an so on to column 4 - point 4
% a third index identifies the Panel.
% And so PanelPoints(i,j,k) contains the i-th coordinate (x,y, or z)
% of the j-th point (FL, FR, AR, or AL) of the k-th Panel.
%
% Max number of spanwise strips and chordwise divisions
% per any panel:

%NSmax=30;    % This was previously inputted here. It is now read from the input.json file.                                 
%NCmax=30;    % This was previously inputted here. It is now read from the input.json file.                                     
%
% Check
%
for ipanel=1:nPanels % Loop over all panels
    ns=data_panel(ipanel,13); % ns - number of divisions spanwise of panel ipanel
    nc=data_panel(ipanel,14); % nc - number of divisions chordwise of panel ipanel

    if ((ns > NSmax) || (nc > NCmax))
        disp(' Too many divisions per panel ');
        return;
    else
    end
end % end of loop over panels



% Check y coordinates (the code can only handle panels with:
% yFR=yAR and yAL=yFL) So check that this is true.

for ipanel=1:nPanels % Loop over all panels
    yFRcheck=data_panel(ipanel,5); 
    yARcheck=data_panel(ipanel,8); 
    yFLcheck=data_panel(ipanel,2); 
    yALcheck=data_panel(ipanel,11); 

    if yFRcheck ~= yARcheck || yFLcheck ~= yALcheck
        disp([' The neccesary condition of yFR=yAR and yAL=yFL is not ' ...
            'satisfied for at least one panel defined by user ']);
        return;
    else
    end
end % end of loop over panels


%Also check that zFR=zAR and zFL=zAL for all panels
%This is necessary to ensure that non-vertical panels are planar
%And for vertical panels, this ensures that the left and right
%sides of the panel are aligned with the x-axis. 

for ipanel=1:nPanels % Loop over all panels
    zFRcheck=data_panel(ipanel,6); 
    zARcheck=data_panel(ipanel,9); 
    zFLcheck=data_panel(ipanel,3); 
    zALcheck=data_panel(ipanel,12); 

    if zFRcheck ~= zARcheck || zFLcheck ~= zALcheck
        disp([' The neccesary condition of zFR=zAR and zAL=zFL is not ' ...
        'satisfied for at least one panel defined by user ']);
        return;
    else
    end
end 

%Check that aspect ratio of each box is less than or equal to 2.
%This is necessary because the 1972 implementation of DLM limits aspect
%ratios.

for ipanel=1:nPanels
    % if max(norm(FR-FL),norm(AR-AL))/min(norm(FL-AL),norm(FR-AR))  
    % This is a conservative estimate of aspect ratio of each box.
    if max(norm(data_panel(ipanel,4:6)-data_panel(ipanel,1:3)),norm(data_panel(ipanel,7:9)-data_panel(ipanel,10:12)))/min(norm(data_panel(ipanel,1:3)-data_panel(ipanel,10:12)),norm(data_panel(ipanel,4:6)-data_panel(ipanel,7:9)))>2*data_panel(ipanel,14)/data_panel(ipanel,13)
        error('Aspect ratios of individual boxes are too large, make mesh finer such that aspect ratio of each box is below 2');
    end
end 

%Check whether you have any co-planar panels.
% Give warning if you do.
% This is because co-planar panels that are one behind anothere must be divided such
%that the mesh boxes in these panels are aligned as to avoid the
%trailing vortex from one panel hitting the collocation point on the other.

for ipanel=1:nPanels
    for jpanel=1:nPanels
        %Check is the two panels are parallel by comparing their normals in
        %the next if statement
        normal_ipanel=cross(data_panel(ipanel,4:6)-data_panel(ipanel,1:3),data_panel(ipanel,4:6)-data_panel(ipanel,7:9));
        normal_ipanel=normal_ipanel/norm(normal_ipanel);
        normal_jpanel=cross(data_panel(jpanel,4:6)-data_panel(jpanel,1:3),data_panel(jpanel,4:6)-data_panel(jpanel,7:9));
        normal_jpanel=normal_jpanel/norm(normal_jpanel);
        if abs(normal_jpanel(1)-normal_ipanel(1))<0.00001 && abs(normal_jpanel(2)-normal_ipanel(2))<0.00001 && abs(normal_jpanel(3)-normal_ipanel(3))<0.00001
            %Now check whether ith panel and jth panel are coplanar for
            %vertical and non-vertical panels 
            if abs(abs(normal_ipanel(2))-1)<0.00001  %Compare abs(y coordinate) of normalized normal vector to 1 to see if panel is vertical 
                if abs(data_panel(ipanel,5)-data_panel(jpanel,5))<0.00001 % compare the y coordinate of FR corner of both vertical panels to see if coplanar
                    warning('This is a warning that there are coplanar panels. If these panels are one behind another in the streamwise direction, do not forget to make sure that the streamwise edges of the mesh boxes in the two panels are aligned!');
                end
            else % Panel is non-vertical
                if abs(data_panel(ipanel,6)-data_panel(jpanel,6))<0.00001 % compare the z coordinate of FR corner of both non-vertical panels to see if coplanar
                    warning('This is a warning that there are coplanar panels. If these panels are one behind another in the streamwise direction, do not forget to make sure that the streamwise edges of the mesh boxes in the two panels are aligned!');
                end
            end
        end
    end
end
        

% Allocate memory for PanelPoints
PanelPoints=zeros(3,4,nPanels);
% Each column contains the x,y,z coordinates of a vertex of the panel 
% with column 1 - point FL up to column 4 with point AL. The 3rd index
% tells us what panel it is.
%
nBoxesTotal=0;  % Total Number of aerodynamic boxes (to be found next)
%
% Number of chordwise divisions for each panel
NchordDivisions=zeros(nPanels,1); % a vector of length nPanels (just initilizing here)
% Number of spanwise divisions for each panel
NspanDivisions=zeros(nPanels,1); % a vector of length nPanels (just initilizing here)
% Number of boxes for each panel
NboxesPerPanel=zeros(nPanels,1); % a vector of length nPanels (just initilizing here)
%
% Panel Definitions
%
for ip=1:nPanels %                   Loop over all panels
    %
    % For each panel enter the x,y,z coords of its four vertex points into
    % the array PanelPoints
    %
    % FL (fwd-left) into column 1
    PanelPoints(1,1,ip)=data_panel(ip,1);
    PanelPoints(2,1,ip)=data_panel(ip,2);
    PanelPoints(3,1,ip)=data_panel(ip,3);
    %
    % FR (fwd-right) into column 2
    PanelPoints(1,2,ip)=data_panel(ip,4);
    PanelPoints(2,2,ip)=data_panel(ip,5);
    PanelPoints(3,2,ip)=data_panel(ip,6);
    %
    % AR (aft-right) into column 3
    PanelPoints(1,3,ip)=data_panel(ip,7);
    PanelPoints(2,3,ip)=data_panel(ip,8);
    PanelPoints(3,3,ip)=data_panel(ip,9);
    %
    % AL (aft-left) into column 4
    PanelPoints(1,4,ip)=data_panel(ip,10);
    PanelPoints(2,4,ip)=data_panel(ip,11);
    PanelPoints(3,4,ip)=data_panel(ip,12);
    %
    % Divisions chordwise:
    NchordDivisions(ip)=data_panel(ip,13); % nc for the panel
    % Divisions spanwise:
    NspanDivisions(ip)=data_panel(ip,14); % ns for the panel
    % number of boxes on panel ip
    NboxesPerPanel(ip)= NchordDivisions(ip)*NspanDivisions(ip); % nc x ns 
    %
    % accumulating total number of boxes: nBoxesTotal
    %
    nBoxesTotal=nBoxesTotal+NboxesPerPanel(ip);
    %
    % End the loop on all panels
end
%
% --------------------------------------------------------
% Now to the generation of aerodynamic boxes
% --------------------------------------------------------
%
% Allocate arrays for aerodynamic boxes and their doublets
%
% Array "Boxes" contains the x,y,z coordinates of corner points of boxes
% plus identification indices
%
% first index 1-3 : location of point (x,y,z)
% second index 1-4: identity of corner point for the box
% index - box number on the panel: idBox
%
% Array iBoxesTable contains the location of box i,j of panel k
% in the list of all boxes for the configuration
%
%
Boxes=zeros(3,4,nBoxesTotal); % Allocate memory for array Boxes
% BoxArea is defined as a column vector
BoxArea=zeros(nBoxesTotal,1); % Allocate memory for vector of box areas
%
ProjectedBoxArea=zeros(nBoxesTotal,1); % a vector
% Cosines of Dihedral angles of Boxes
CosDihedral=zeros(nBoxesTotal,1); % a vector
% unit normals to boxes 
% boxes that are vertical are consistent with the way normals are
% defined for non vertical boxes. That is, for verrtical boxes
% the normal is in the -y direction
normals=zeros(3,nBoxesTotal); 
% The normal to each box is stored in the column for that box in array normals
% 
% More storage allocation for arrays
%
iBoxesTable=zeros(NSmax,NCmax,nPanels); % nc and nc for each panel
iBoxBelongsToPanel=zeros(nBoxesTotal); % for each box - what panel it belongs to
iBoxOnPanelChordwise=zeros(nBoxesTotal); % for each box - its place chordwise
iBoxOnPanelSpanwise=zeros(nBoxesTotal); % for each box - its place spanwise
%
% Coefficients of LE and TE equations for panels (see later
% xLE0, xLE1, etc.)
%
% Each panel has equations for its LE and TE in the form:
%
% Prepare array storage for the coefficients of the LE and TE lines.
%
LEx0=zeros(nPanels,1);
LEx1=zeros(nPanels,1);
LEz0=zeros(nPanels,1);
LEz1=zeros(nPanels,1);
TEx0=zeros(nPanels,1);
TEx1=zeros(nPanels,1);
TEz0=zeros(nPanels,1);
TEz1=zeros(nPanels,1);
%
%
% Find x,y,z coordinates of corner points for all boxes
%
% Loop over panels (trapezoids)
%
% Aerodynamic Box identification index: idBox
%
idBox=0; % Box identifier
%
for ipanel=1:nPanels % Loop over the panel
    %
    % For each panel - create the LE and TE equations
    %
    % Equation defining the x,y,z coordinates of points on the leading edge
    % of panel ipanel - straight line from point 1 (FL) to point 2 (FR)
    % L is for Left. R is for Right
    %
    %
    % The leading edge (LE) based on FL and FR
    %
    yL=PanelPoints(2,1,ipanel);
    yR=PanelPoints(2,2,ipanel);
    xL=PanelPoints(1,1,ipanel);
    xR=PanelPoints(1,2,ipanel);
    zL=PanelPoints(3,1,ipanel);
    zR=PanelPoints(3,2,ipanel);

    %
    % (x-xL)/(xR-xL)=(y-yL)/(yR-yL) % The line connecting the L and R
    % points:
    %
    % The x=x(y) LE equation:
    % x=xL+(xR-xL)/(yR-yL)*(y-yL)=xLE0+xLE1*y
    %
    % -------------------------------------------------
    % Non vertical Panels
    % -------------------------------------------------
    %
    if abs(yR-yL)>1.d-5 
    % Check that the R and L points do not have the same y location.
    
    % Check that user has defined the corner points such that yR>yL:
    % This is necessary to avoid getting negative distances in deltaY=(yR-yL)/(NS);
        if yL>yR
            disp(['The corner points of all non-vertical panels must be ' ...
                'defined such that yR>yL (the direction "right" is defined as ' ...
                'the direction with the more positive y value)'])
            return;
        else
        end

    % The x=x(y) LE equation:
        xLE1=(xR-xL)/(yR-yL);
        xLE0=xL-xLE1*yL;
        % store in array:
        LEx0(ipanel)=xLE0;
        LEx1(ipanel)=xLE1;
        %
        % Similarly for z:
        % The z=z(y) equation for the LE:
        % z=zL+(zR-zL)/(yR-yL)*(y-yL)=zLE0+zLE1*y
        %
        zLE1=(zR-zL)/(yR-yL);
        zLE0=zL-zLE1*yL;
        %
        % store in array:
        LEz0(ipanel)=zLE0;
        LEz1(ipanel)=zLE1;
        %
        %
        % Now to the TE equations x=x(y) and z=z(y):
        %
        % Equation defining the x,y,z coordinates of points on the trailing edge
        % of panel ipanel - straight line from point 4 (AL) to point 3 (AR)
        %
        yL=PanelPoints(2,4,ipanel);
        yR=PanelPoints(2,3,ipanel);
        xL=PanelPoints(1,4,ipanel);
        xR=PanelPoints(1,3,ipanel);
        zL=PanelPoints(3,4,ipanel);
        zR=PanelPoints(3,3,ipanel);
        %
        % (x-xL)/(xR-xL)=(y-yL)/(yR-yL)
        %
        % The x=x(y) equation for the TE:
        %
        % x=xL+(xR-xL)/(yR-yL)*(y-yL)=xTE0+xTE1*y
        %
        xTE1=(xR-xL)/(yR-yL);
        xTE0=xL-xTE1*yL;
        %
        % store in array:
        TEx0(ipanel)=xTE0;
        TEx1(ipanel)=xTE1;
        %
        % Similarly for z:
        %
        % The z=z(y) equation for the TE:
        %
        % z=zL+(zR-zL)/(yR-yL)*(y-yL)=zTE0+zTE1*y
        %
        zTE1=(zR-zL)/(yR-yL);
        zTE0=zL-zTE1*yL;
        %
        % store in array:
        TEz0(ipanel)=zTE0;
        TEz1(ipanel)=zTE1;
        %
        % We now have equations for x and z as functions of y for the LE
        % and TE.
        %
        % Loop over spanwise strips of boxes for this panel
        %
        NS=NspanDivisions(ipanel);
        NC=NchordDivisions(ipanel);
        %
        % Width of each strip
        %
        deltaY=(yR-yL)/(NS);
        %
        % left y and right y coordinates of first spanwise strip
        yleft=yL-deltaY;   % just used for initialization
        yright=yL;         % just used for initialization
        %
        % x coordinate of the left forward (LF) point
        %xstart=PanelPoints(1,1,ipanel); -Line commented out because
        %variable was defined but never used 
        %
        for ispanwise=1:NS % Loop over the spanwise strips
            %
            % strip ispanwise
            %
   % move from strip to strip by augmenting their left and right locations
            yleft=yleft+deltaY; 
            yright=yright+deltaY;
            %
            % x and z coordinates of LE and TE points on the left line from
            % LE to TE
            %
            xLEleft=xLE0+xLE1*yleft;
            xTEleft=xTE0+xTE1*yleft;
            zLEleft=zLE0+zLE1*yleft;
            zTEleft=zTE0+zTE1*yleft;
            %
            % x and z coordinates of LE and TE points on the right line
            % from LE to TE
            %
            xLEright=xLE0+xLE1*yright;
            xTEright=xTE0+xTE1*yright;
            zLEright=zLE0+zLE1*yright;
            zTEright=zTE0+zTE1*yright;
            %
            % We now have the x,z locations of the points on the LE and TE
            % along the two streamwise lines that define the left and right
            % sides of the strip (strip ispanwise)
            %
            % We proceed to divide the strip into boxes in the chordwise
            % (x) direction.
            %
            %
            %
            % x and z step size chordwise on left:
            %
            dXleft=(xTEleft-xLEleft)/NC;
            dZleft=(zTEleft-zLEleft)/NC;
            %
            % step size chordwise on right:
            %
            dXright=(xTEright-xLEright)/NC;
            dZright=(zTEright-zLEright)/NC;
            %
            % Loop over chordwise boxes for this spanwise strip
            %
            xLEleft=xLEleft-dXleft;
            xLEright=xLEright-dXright;
            zLEleft=zLEleft-dZleft;
            zLEright=zLEright-dZright;
            %
            for ichordwise=1:NC 
                %
                % Count the box number in the total list of all boxes
                %
                idBox=idBox+1;
                %
                xLEleft=xLEleft+dXleft;
                xLEright=xLEright+dXright;
                zLEleft=zLEleft+dZleft;
                zLEright=zLEright+dZright;
                %
                % Array Boxes contains the x,y,z coordinates of
                % corner points of boxes
                % plus identification indices
                %
                % first index: 1 or 2 or 3 : location of point (x,y,z)
                % second index: 1 or 2 or 3 or 4: identity of corner point for the box
                % third index - box number idBox in the set of all boxes
                %
                % Array iBoxesTable contains the location of box i,j of panel k
                % in the list of all boxes for the configuration
                % As a reminder of how these arrays are structured:
                % Boxes=zeros(3,4,nBoxMax)             
                % iBoxesTable=zeros(NSmax,NSmax,nPanels)
                %
                iBoxesTable(ispanwise,ichordwise,ipanel)=idBox;
                %
                % Store, for box idBox, the identity of the panel
                % and, on that panel, its chordwise and spanwise positions
                %
                 iBoxBelongsToPanel(idBox)=ipanel;
                 iBoxOnPanelChordwise(idBox)=ichordwise;
                 iBoxOnPanelSpanwise(idBox)=ispanwise;
                %
                x1=xLEleft;
                y1=yleft;
                z1=zLEleft;
                %
                x2=xLEright;
                y2=yright;
                z2=zLEright;
                %
                x3=xLEright+dXright;
                y3=yright;
                z3=zLEright+dZright;
                %
                x4=xLEleft+dXleft;
                y4=yleft;
                z4=zLEleft+dZleft;
                %
                Boxes(1,1,idBox)=x1;
                Boxes(2,1,idBox)=y1;
                Boxes(3,1,idBox)=z1;
                Boxes(1,2,idBox)=x2;
                Boxes(2,2,idBox)=y2;
                Boxes(3,2,idBox)=z2;
                Boxes(1,3,idBox)=x3;
                Boxes(2,3,idBox)=y3;
                Boxes(3,3,idBox)=z3;
                Boxes(1,4,idBox)=x4;
                Boxes(2,4,idBox)=y4;
                Boxes(3,4,idBox)=z4;
   %
   % Account for dihedral of box
   %
                L=sqrt((y2-y1)^2+(z2-z1)^2);
                CosDihedral(idBox)=(y2-y1)/L;
                BoxArea(idBox)=...
                abs((y2-y1)/2.*(x3-x2+x4-x1)/CosDihedral(idBox));
                ProjectedBoxArea(idBox)=...
                abs((y2-y1)/2.*(x3-x2+x4-x1)); 
    %
    % Find unit normal to box area 
    % 
                for ii=1:3
                    R1(ii)=Boxes(ii,4,idBox)-Boxes(ii,1,idBox);%vector from point 1 to pt 4 (along root)
                    R2(ii)=Boxes(ii,2,idBox)-Boxes(ii,1,idBox);%vector from pt 1 to pt2 (along LE) 
                end
                % Instead of the loop over vector components 1,2,3 above, this can be
                % written more concisely as:
                %R1(:)=Boxes(:,4,idBox)-Boxes(:,1,idBox);%vector from point 1 to pt 4 (along root)
                %R2(:)=Boxes(:,2,idBox)-Boxes(:,1,idBox);%vector from pt 1 to pt2 (along LE) 
                %
                R1xR2=cross(R1,R2); % vector product of R1 x R2 = vector perpendicular to the box
                Aux=R1xR2/norm(R1xR2,2); % The length of that vector
                % Generating a unit vector that is perpendicular to the surface, pointing
                % up
                normals(1,idBox)=Aux(1);
                normals(2,idBox)=Aux(2);
                normals(3,idBox)=Aux(3);
                % We can write this more concisely: normals(:,idBox)=Aux(:)
                %
                % end loop over chordwise boxes
                %
            end
            %
            % end loop over spanwise divisions
            %
        end
        % End branch of if statement that checked if a panel is vrertical.
        % In the branch above the panel is not vertical.
        % after the else command, vertical panels are treated.
    else
    %
    % --------------------------------------------------------------
    % The case of a Vetical Panel (No Equivalent Dihedral Allowed)
    % --------------------------------------------------------------

    % The LE and TE lines for x are now functions of z
    % For the LE:
    % x=xLE0+xLE1*z
    %
    % xL,xR, etc. for the LE were already created at the top of the loop on panels
    % (ipanel) before the check (if statement) on whether the panel is
    % vertical was done.
    %
        xLE1=(xR-xL)/(zR-zL);
        xLE0=xL-xLE1*zL;
        % store in array:
        LEx0(ipanel)=xLE0;
        LEx1(ipanel)=xLE1;
        %
        % Equation defining the x,y,z coordinates of points on the trailing edge
        % of panel ipanel - straight line from point 4 (AL) to point 3 (AR)
        %
        yL=PanelPoints(2,4,ipanel);
        yR=PanelPoints(2,3,ipanel);
        xL=PanelPoints(1,4,ipanel);
        xR=PanelPoints(1,3,ipanel);
        zL=PanelPoints(3,4,ipanel);
        zR=PanelPoints(3,3,ipanel);

        if zL>zR  % Check that zR>zL (otherwise the code will give 
            % negative distances while computing deltaZ=(zR-zL)/(NS); )
            disp(['The corner points of all vertical panels must be ' ...
                'defined such that zR>zL (the direction "right" is ' ...
                'defined as the direction with the more positive z value)'])
            return;
        else
        end


        % z would replace y in the following equations:
        %
        % (x-xL)/(xR-xL)=(y-yL)/(yR-yL)
        %
        % x=xL+(xR-xL)/(yR-yL)*(y-yL)=xTE0+xTE1*y
        %
        xTE1=(xR-xL)/(zR-zL);
        xTE0=xL-xTE1*zL;
        % x=xTE0+xTE1*z
        %
        % store in array:
        TEx0(ipanel)=xTE0;
        TEx1(ipanel)=xTE1;
        %
        % Loop over spanwise strips of boxes for this panel
        %
        NS=NspanDivisions(ipanel);
        NC=NchordDivisions(ipanel);
        %
        % Width of each strip
        %
        deltaZ=(zR-zL)/(NS);
        %
        % left y and right y coordinates of first spanwise strip
        zleft=zL-deltaZ;
        zright=zL;
        %
        % x coordinate of the left forward (LF) point
        % xstart=PanelPoints(1,1,ipanel); Line commented out because
        % variable was defined but never used. 
        %
        %
        for ispanwise=1:NS % Loop over spanwise strips
            %
            % strip ispanwise
            %
            % Augment the left and right sides of the strip (bottom and
            % top sides) by deltaZ
            %
            zleft=zleft+deltaZ;
            zright=zright+deltaZ;
            %
            % x coordinates of LE and TE points on the left line
            % of the strip
            xLEleft=xLE0+xLE1*zleft;
            xTEleft=xTE0+xTE1*zleft;
            %
            % x coordinates of LE and TE points on the right line
            % for the strip:
            %
            xLEright=xLE0+xLE1*zright;
            xTEright=xTE0+xTE1*zright;
            % x step size chordwise on left:
            %
            dXleft=(xTEleft-xLEleft)/NC;
            %
            % step size chordwise on right:
            %
            dXright=(xTEright-xLEright)/NC;
            %
            % Loop over chordwise boxes for this spanwise strip
            %
            xLEleft=xLEleft-dXleft;
            xLEright=xLEright-dXright;
            %
            for ichordwise=1:NC % Loop over the chordwise boxes
                %
                % Count the box number in the total list of boxes
                %
                idBox=idBox+1; % Keep counting boxes for the full set of boxes
                %
                xLEleft=xLEleft+dXleft;
                xLEright=xLEright+dXright;
                %
                % Array Boxes contains the x,y,z coordinates of
                % corner points of boxes
                % plus identification indices
                %
                % first index 1 or 2 or 3 : location of point (x,y,z)
                % second index 1 or 2 or 3 or 4: identity of corner point for the box
                % third index - box number in the list of all boxes: idBox
                % Boxes=zeros(3,4,nBoxMax)
                %
                % Array iBoxesTable contains the location of box i,j of panel k
                % in the list of all boxes for the configuration
                % iBoxesTable=zeros(NSmax,NSmax,nPanels)
                %
                iBoxesTable(ispanwise,ichordwise,ipanel)=idBox;
                %
                % Note: All vertical boxes are at yL=yR. 
                % It is important to identify the side: the direction of 
                % the normal to the box.
                % The normal is in the -y direction (see later)
                %
                x1=xLEleft;
                y1=yL;
                z1=zleft;
                %
                x2=xLEright;
                y2=yL;
                z2=zright;
                %
                x3=xLEright+dXright;
                y3=yL;
                z3=zright;
                %
                x4=xLEleft+dXleft;
                y4=yL;
                z4=zleft;
                %
                Boxes(1,1,idBox)=x1;
                Boxes(2,1,idBox)=y1;
                Boxes(3,1,idBox)=z1;
                Boxes(1,2,idBox)=x2;
                Boxes(2,2,idBox)=y2;
                Boxes(3,2,idBox)=z2;
                Boxes(1,3,idBox)=x3;
                Boxes(2,3,idBox)=y3;
                Boxes(3,3,idBox)=z3;
                Boxes(1,4,idBox)=x4;
                Boxes(2,4,idBox)=y4;
                Boxes(3,4,idBox)=z4;
                %
                % Only vertical boxes allowed.
                % That is, the dihedral angle is 90 deg.
                CosDihedral(idBox)=0.0;
                normals(1,idBox)=0.;
                normals(2,idBox)=-1.;
                normals(3,idBox)=0.;
                % Note: The normal points in the -y direction
                % Note: It is as if we defined the vertical panel while
                % lying in the xy plane and then rotated it 90 degrees
                % to a vertical location
                ProjectedBoxArea(idBox)=0.0; % projection on the xy plane
                BoxArea(idBox)=abs((z2-z1)/2.*(x3-x2+x4-x1));
                %
                % end loop over chordwise boxes
                %
            end
            %
            % end loop over spanwise divisions
            %
        end
        %
        % end if that checked vertical position or not
    end
    %
    % end loop over all panels to define geometry of all boxes
    %
end
%
nBoxes=idBox; % The total number of boxes.
%
% ---------------------------------------------------------------
% At this stage all aerodynamic boxes have been defined and the following
% arrays have been generated:
%
% For the panels:
% data_panel(nPanels,14)
% PanelPoints(3,4,nPanels)
% NchordDivisions(nPanels)
% NboxesPerPanel(nPanels)
% iBoxesTable(NSmax,NCmax,nPanels) - nc and nc for each panel
%
% For the boxes:
% Boxes(3,4,nBoxesMax)
% ProjectedBoxArea(nBoxesMax)
% BoxArea(nBoxesMax)
% CosDihedral(nBoxesMax)
% iBoxBelongsToPanel=zeros(nBoxesMax); % for each box - what panel it belongs to
% iBoxOnPanelChordwise=zeros(nBoxesMax); % for each box - its place chordwise
% iBoxOnPanelSpanwise=zeros(nBoxesMax); % for each box - its place spanwise
%
%
% ---------------------------------------------------------------
%
%
% ----------------------------------------------------------------
% Plotting the mesh
% ----------------------------------------------------------------
%
% Plot vertex points for all boxes
%
if(plotbox==1) % check if the control index for plotting boxes is 1. 
    %
    for ipanel=1:nPanels % Loop over all panels
        %
        NS=NspanDivisions(ipanel);
        NC=NchordDivisions(ipanel);
        % yL and yR coordinates
        yL=PanelPoints(2,4,ipanel);
        yR=PanelPoints(2,3,ipanel);
        % aft left
        xAL=PanelPoints(1,4,ipanel);
        zAL=PanelPoints(3,4,ipanel);
        % aft right
        xAR=PanelPoints(1,3,ipanel);
        zAR=PanelPoints(3,3,ipanel);
        %
        % fwd left
        xFL=PanelPoints(1,1,ipanel);
        zFL=PanelPoints(3,1,ipanel);
        % fwd right
        xFR=PanelPoints(1,2,ipanel);
        zFR=PanelPoints(3,2,ipanel);
        %
        % Division into boxes at the root and the tip
        %
        dxRoot=(xAL-xFL)/NC;
        dzRoot=(zAL-zFL)/NC;
        dxTip=(xAR-xFR)/NC;
        dzTip=(zAR-zFR)/NC;
        %
        %
        % First: Non vertical panels
        %
        if abs(yR-yL)> 1e-7 
            %
            deltaY=(yR-yL)/NS;
            %
            % loop over spanwise lines going, each, from yL to yR
            %
            for iline=1:(NC+1)
                %
                % left point on line
                yy(1)=yL;
                xx(1)=xFL+(iline-1)*dxRoot;
                zz(1)=zFL+(iline-1)*dzRoot;
                % right point on line
                yy(2)=yR;
                xx(2)=xFR+(iline-1)*dxTip;
                zz(2)=zFR+(iline-1)*dzTip;
                %
                plot3(xx,yy,zz,'-'); %                 Plot3 line by line
                hold on;
            end   % end loop on left-to-right lines
            axis('equal')
            %
            clear xx;
            clear yy;
            clear zz;
            %
            % loop over chordwise lines
            %
            for iline=1:(NS+1)
                %
                % fwd point on line ("sets" on LE of panel)
                yy(1)=yL+(iline-1)*deltaY;
                xx(1)=LEx0(ipanel)+LEx1(ipanel)*yy(1);
                zz(1)=LEz0(ipanel)+LEz1(ipanel)*yy(1);
                % aft point on line ("seats" on TE of panel)
                yy(2)=yy(1);
                xx(2)=TEx0(ipanel)+TEx1(ipanel)*yy(1);
                zz(2)=TEz0(ipanel)+TEz1(ipanel)*yy(1);
                %
                plot3(xx,yy,zz,'-'); %               Plot 3 chordwise lines
            end
        %    
        else
            % ---------------------------
            % The Case of Vertical Panels
            % ---------------------------
            %
            deltaZ=(zFR-zFL)/NS;
            %
            for iline=1:(NC+1) % spanwise lines at different z distances
                %
                % "left" point on line
                yy(1)=yL;  % Note that in the case of vertical panels yL=yR
                xx(1)=xFL+(iline-1)*dxRoot;
                zz(1)=zFL+(iline-1)*dzRoot;
                % "right" point on line
                yy(2)=yL;
                xx(2)=xFR+(iline-1)*dxTip;
                zz(2)=zFR+(iline-1)*dzTip;
                %
                plot3(xx,yy,zz,'-'); %            Plot3
                hold on;
            end             % end loop on chordwise lines
            axis('equal')
            %
            clear xx;
            clear yy;
            clear zz;
            %
            % loop over chordwise lines
            %
            for iline=1:(NS+1)
                %
                % fwd point on line
                yy(1)=yL;
                zz(1)=zFL+(iline-1)*deltaZ;
                xx(1)=LEx0(ipanel)+LEx1(ipanel)*zz(1);
                % aft point on line
                yy(2)=yL;
                zz(2)=zz(1);
                xx(2)=TEx0(ipanel)+TEx1(ipanel)*zz(2);
                %
                plot3(xx,yy,zz,'-');
            end     % end loop on spanwise lines
        end % end of if that checks for verical panels
    end     % end of loop on panels
end % end of if that checks whether to plot or not
%
% -----------------------------------------------
% Beginning of the Doublet Lattice Process
% -----------------------------------------------
%
% Loop over all aerodynamic boxes
% to generate force and downwash points
% as well as 1/4 chord left and right tip points
%
BoxSendPoint1=zeros(nBoxesTotal,3);
BoxSendPoint2=zeros(nBoxesTotal,3);
BoxSendPoint3=zeros(nBoxesTotal,3);
BoxPointDownwash=zeros(nBoxesTotal,3);
%
% The deltaY width of the boxes stored in a vector
%
BoxWidth=zeros(nBoxesTotal,1);
BoxPointChord=zeros(nBoxesTotal,1);
%
for idBox=1:nBoxesTotal % Loop over all boxes
    %
    % create for each box the coordinates of points 1 & 2 (the
    % corners of the doublet), the coordinates of the sending point,
    % and the coordinates of the receiving point
    %
    % Find the 4 vertex points of the box
    %
    x1=Boxes(1,1,idBox);
    y1=Boxes(2,1,idBox);
    z1=Boxes(3,1,idBox);
    x2=Boxes(1,2,idBox);
    y2=Boxes(2,2,idBox);
    z2=Boxes(3,2,idBox);
    x3=Boxes(1,3,idBox);
    y3=Boxes(2,3,idBox);
    z3=Boxes(3,3,idBox);
    x4=Boxes(1,4,idBox);
    y4=Boxes(2,4,idBox);
    z4=Boxes(3,4,idBox);
    %
    % Distinguish between non-vertical and vertical boxes
    %
    if abs(y2-y1)>1e-7  % Check if the box is vertical or not
        % The non-vertical box case
        % Should we correct for Dihedral Effect here ?????????????????????????
        BoxWidth(idBox)=y2-y1;
    else
        % The vertical box case
        BoxWidth(idBox)=z2-z1;
    end
    %
    % The sending points is at 0.25 chord inboard[Point1], outboard[Point2]
    % and center[Point3] load point
    % two points - at the left and right of the 1/4 chord line
    % and one point (the load point) at the center of the 1/4 chord line.
    %
    if abs(y2-y1)>1e-7  % Check if the box is vertical or not
        %
        % The two ends of the line of the vortex
        %
        BoxSendPoint1(idBox,1)=x1+0.25*(x4-x1);
        BoxSendPoint1(idBox,2)=y1;
        BoxSendPoint1(idBox,3)=z1+0.25*(z4-z1);
        %
        BoxSendPoint2(idBox,1)=x2+0.25*(x3-x2);
        BoxSendPoint2(idBox,2)=y2;
        BoxSendPoint2(idBox,3)=z2+0.25*(z3-z2);
        %
        % The location of the force point (midway between Point1 and Point2)
        for ii=1:3
                BoxSendPoint3(idBox,ii)=(BoxSendPoint1(idBox,ii)...
                                        +BoxSendPoint2(idBox,ii))/2.;
        end
        %
        % The downwash point is at 0.75 chord between (spanwise) points located
        % at the 3/4 chord of the left and right boundaries of the box
        %
        BoxPointDownwash(idBox,1)=0.5*(x1+0.75*(x4-x1)+x2+0.75*(x3-x2));
        BoxPointDownwash(idBox,2)=0.5*(y1+y2);
        BoxPointDownwash(idBox,3)=0.5*(z1+0.75*(z4-z1)+z2+0.75*(z3-z2));
        %
    else
        %
        % Vertical boxes
        %
        % The two ends of the line of the vortex
        %
        BoxSendPoint1(idBox,1)=x1+0.25*(x4-x1);
        BoxSendPoint1(idBox,2)=y1;
        BoxSendPoint1(idBox,3)=z1;
        %
        BoxSendPoint2(idBox,1)=x2+0.25*(x3-x2);
        BoxSendPoint2(idBox,2)=y2;
        BoxSendPoint2(idBox,3)=z2;
        %
        % The location of the force point (midway between Point1 and Point2)
        for ii=1:3
                BoxSendPoint3(idBox,ii)=(BoxSendPoint1(idBox,ii)...
                                        +BoxSendPoint2(idBox,ii))/2.;
        end
        %
        % The downwash point is at 0.75 chord between (spanwise) points located
        % at the 3/4 chord of the left and right boundaries of the box
        %
        BoxPointDownwash(idBox,1)=0.5*(x1+0.75*(x4-x1)+x2+0.75*(x3-x2));
        BoxPointDownwash(idBox,2)=0.5*(y1+y2);
        BoxPointDownwash(idBox,3)=0.5*(z1+z2);
        %
        % end of if checking if boxes are vertical
        %
    end
    %
    % chord of box
    %
    BoxPointChord(idBox)=((x4-x1)+(x3-x2))/2;
%
end
%
% End of loop over all boxes to determine their 1/4 chord corner points
% and their load and downwash points
%
% We now have the coordinates of:
% Points 1&2 of each box: The ends of the vortex line at 0.25 chord
% BoxSendPoint1(nBoxesTotal,3); 
% BoxSendPoint2(nBoxesTotal,3);
% The force points (where the vorticies are)
% BoxSendPoint3(nBoxesTotal,3);
% The downwash points, where normalwash BC are enforced
% BoxPointDownwash(nBoxesTotal,3);
%
%
% --------------------------------------------
% ---------------------------
% More plotting
% ---------------------------
% A reminder:
% Control of plotting (zeros or 1s):
% element 1 - plot panels & Aero Boxes
% element 2 - plot vortex edge points on boxes
% Element 3 - plot force points and downwash on boxes
%
%
% plotbox=plot_control(1);
% plotVortex=plot_control(2);
% plotforceanddownwash=plot_control(3);
%
% plot load and downwash points.
%
    for ii=1:nBoxesTotal % Loop over all boxes
        %
        xx1(ii)=BoxSendPoint1(ii,1);
        yy1(ii)=BoxSendPoint1(ii,2);
        zz1(ii)=BoxSendPoint1(ii,3);
        %
        xx2(ii)=BoxSendPoint2(ii,1);
        yy2(ii)=BoxSendPoint2(ii,2);
        zz2(ii)=BoxSendPoint2(ii,3);
        %
        xxForce(ii)=BoxSendPoint3(ii,1);
        yyForce(ii)=BoxSendPoint3(ii,2);
        zzForce(ii)=BoxSendPoint3(ii,3);
        %
        xxD(ii)=BoxPointDownwash(ii,1);
        yyD(ii)=BoxPointDownwash(ii,2);
        zzD(ii)=BoxPointDownwash(ii,3);
    end
    %
    %
 if (plotforceanddownwash==1)
 hold on
 plot3(xxForce,yyForce,zzForce,'gd',...
       xxD,yyD,zzD,'ro');
    axis('equal')
 
else
 end
 if (plotVortex==1)
 hold on
 plot3(xx1,yy1,zz1,'gd',xx2,yy2,zz2,'ro');
    axis('equal')
else
 end






%%%%%THIS PART TRANSFERS THE GEOMETRY INTO A .JSON SO THAT IT CAN BE READ BY THE PYTHON
%%%%%PROGRAM OF DLR.

%%%%%BELOW IS USED FOR RUNNING THE PYTHON CODE.

%%%%To run the Python code, do this:

% Step 1: Create the struct (dictionary-like object)
aerogrid = struct();
aerogrid.n = length(BoxPointChord);
aerogrid.l = BoxPointChord;
aerogrid.N = normals';
aerogrid.offset_j = BoxPointDownwash';
aerogrid.offset_l = BoxSendPoint3';
aerogrid.offset_P1 = BoxSendPoint1';
aerogrid.offset_P3 = BoxSendPoint2';

aerogrid.A = BoxArea;


% Step 2: Convert the struct to a JSON string
jsonStr = jsonencode(aerogrid);

% Step 3: Save the JSON string to a .json file
fileID = fopen('aerogrid.json', 'w');  % Open a file for writing
if fileID == -1
    error('Failed to open file for writing.');
end
fwrite(fileID, jsonStr, 'char');  % Write the JSON string to the file
fclose(fileID);  % Close the file