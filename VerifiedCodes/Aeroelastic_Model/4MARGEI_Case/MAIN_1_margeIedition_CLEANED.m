%%This code outputs Mss, Msc, Kss, Css, Q_lst, Q_gust_lst

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
RedFreq_lst=data.RedFreq_lst;  %%defined as wb/U
gust_modeshape=data.gust_modeshape;        %%Use one radian of gust vane deflection at down wash points to define the modeshape
b=data.ref_semichord;                      %Gust modeshape input is a column vector with dimensions: # of aero boxes on gust vanes only x 1
NC=data.n_controlmodes;
n_modes=data.n_modes;
%%Use one radian of surfce deflection at down wash points to define the modeshape

%%Add paths
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'DLM_helperfunctions'));
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'Spline_helperfunctions'));

%%%Get some data from the mesher:
[ProjectedBoxArea,BoxArea,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection] = DLM_Mesh_Updated();


%% The alternative of Above For Marge I Data (Not needed for Marge II)
load('MARGEI_Structural_Modes.mat');
PHI=phi*0.0254;
num_of_nodes=72;
data2=load('MARGEI_nodal_coordinates.mat');
nodal_coordinates=data2.XYZ*0.0254;
b=b*0.0254; %inches to meters conversion of reference chord.                               %<-----------------being divided by 2 at the moment.

ProjectedBoxArea=ProjectedBoxArea*0.0254^2;
BoxArea=BoxArea*0.0254^2;
BoxPointChord=BoxPointChord*0.0254;
BoxPointDownwash=BoxPointDownwash*0.0254;
BoxSendPoint1=BoxSendPoint1*0.0254;
BoxSendPoint2=BoxSendPoint2*0.0254;
BoxSendPoint3=BoxSendPoint3*0.0254;
AD=[2:16,24:31,35:41,48:50]; % AD cannot include point 42,46,47 (as that is part of E), and 51-53 only belong to BC; 
BC=[16:23,31:34,50:53];  
E=[42:43,46:47];       
F=[44:45,48:49];
G=[54:68];                         
H=[65:72];

ADaero=[1:200,401:600];
BCaero=[201:400];
Eaero=[601:700];
Faero=[701:800];
Gaero=[801:900];
Haero=[901:1000];

%% Manipulate the Reduced Frequency List

%!!!!!!!!!!!!
% Note: Reduced Frequency list in user input is omega*b/U but the
% reduced frequency is omega/U in DLM equations,
% not omega*b/U
% So divide RedFreq/b for the input of the DLM_Find_AICs function 
% !!!!!!!!!!!

%If reduced frequency in the user input does not include 0, add it in.
if ~ismember(0, RedFreq_lst)
    RedFreq_lst = [RedFreq_lst; 0]; 
end

RedFreq_lst=sort(RedFreq_lst); 


%% Calculate Aerodynamic Forces except Gust
%%%In this part, we add forces due to all modes (for each frequency)
Q_lst=zeros(n_modes+NC,n_modes+NC,length(RedFreq_lst));    %%list for different frequencies

for ii=1:length(RedFreq_lst)
    RedFreq=RedFreq_lst(ii);
    %Notice division by 'RedFreq/b' in the input of the function below
    %instead of RedFreq.
    [AIC] = DLM_Find_AICs(Mach, RedFreq/b,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection);  
    Q=zeros(n_modes+NC);

    for ishape=1:n_modes+NC
        T_hFE_to_fp=zeros(1600,72);
        T_hFE_to_fp(ADaero,AD)=T_Find(nodal_coordinates(AD,1),nodal_coordinates(AD,2),PHI(AD,ishape),BoxSendPoint3(ADaero,1),BoxSendPoint3(ADaero,2)); % AD spline
        T_hFE_to_fp(BCaero,BC)=T_Find(nodal_coordinates(BC,1),nodal_coordinates(BC,2),PHI(BC,ishape),BoxSendPoint3(BCaero,1),BoxSendPoint3(BCaero,2)); % BC spline
        T_hFE_to_fp(Eaero,E)=T_Find(nodal_coordinates(E,1),nodal_coordinates(E,2),PHI(E,ishape),BoxSendPoint3(Eaero,1),BoxSendPoint3(Eaero,2)); % E spline
        T_hFE_to_fp(Faero,F)=T_Find(nodal_coordinates(F,1),nodal_coordinates(F,2),PHI(F,ishape),BoxSendPoint3(Faero,1),BoxSendPoint3(Faero,2)); % F spline
        T_hFE_to_fp(Gaero,G)=T_Find(nodal_coordinates(G,1),nodal_coordinates(G,2),PHI(G,ishape),BoxSendPoint3(Gaero,1),BoxSendPoint3(Gaero,2)); % G spline
        T_hFE_to_fp(Haero,H)=T_Find(nodal_coordinates(H,1),nodal_coordinates(H,2),PHI(H,ishape),BoxSendPoint3(Haero,1),BoxSendPoint3(Haero,2)); % H spline

        %%Self Induced and Control Forces (not including gust vanes)
        for jshape=1:n_modes+NC
            T_hFE_to_dw=zeros(1600,72);
            T_hFE_to_dw(ADaero,AD)=T_Find(nodal_coordinates(AD,1),nodal_coordinates(AD,2),PHI(AD,jshape),BoxPointDownwash(ADaero,1),BoxPointDownwash(ADaero,2)); % AD spline
            T_hFE_to_dw(BCaero,BC)=T_Find(nodal_coordinates(BC,1),nodal_coordinates(BC,2),PHI(BC,jshape),BoxPointDownwash(BCaero,1),BoxPointDownwash(BCaero,2)); % BC spline
            T_hFE_to_dw(Eaero,E)=T_Find(nodal_coordinates(E,1),nodal_coordinates(E,2),PHI(E,jshape),BoxPointDownwash(Eaero,1),BoxPointDownwash(Eaero,2)); % E spline
            T_hFE_to_dw(Faero,F)=T_Find(nodal_coordinates(F,1),nodal_coordinates(F,2),PHI(F,jshape),BoxPointDownwash(Faero,1),BoxPointDownwash(Faero,2)); % F spline
            T_hFE_to_dw(Gaero,G)=T_Find(nodal_coordinates(G,1),nodal_coordinates(G,2),PHI(G,jshape),BoxPointDownwash(Gaero,1),BoxPointDownwash(Gaero,2)); % G spline
            T_hFE_to_dw(Haero,H)=T_Find(nodal_coordinates(H,1),nodal_coordinates(H,2),PHI(H,jshape),BoxPointDownwash(Haero,1),BoxPointDownwash(Haero,2)); % H spline
            zout_dw=T_hFE_to_dw*PHI(:,jshape);
            
            number_of_aero_panels=length(zout_dw); %%<-- You will need this in the gust force calculation section
           

            zout_div_dx_dw=zeros(1600,1);
            zout_div_dx_dw(ADaero,1)=dTdx_Stickmodel(nodal_coordinates(AD,1),nodal_coordinates(AD,2),PHI(AD,jshape),BoxPointDownwash(ADaero,1),BoxPointDownwash(ADaero,2)); % AD spline
            zout_div_dx_dw(BCaero,1)=dTdx_Stickmodel(nodal_coordinates(BC,1),nodal_coordinates(BC,2),PHI(BC,jshape),BoxPointDownwash(BCaero,1),BoxPointDownwash(BCaero,2)); % BC spline
            zout_div_dx_dw(Eaero,1)=dTdx_Stickmodel(nodal_coordinates(E,1),nodal_coordinates(E,2),PHI(E,jshape),BoxPointDownwash(Eaero,1),BoxPointDownwash(Eaero,2)); % E spline
            zout_div_dx_dw(Faero,1)=dTdx_Stickmodel(nodal_coordinates(F,1),nodal_coordinates(F,2),PHI(F,jshape),BoxPointDownwash(Faero,1),BoxPointDownwash(Faero,2)); % F spline
            zout_div_dx_dw(Gaero,1)=dTdx_Stickmodel(nodal_coordinates(G,1),nodal_coordinates(G,2),PHI(G,jshape),BoxPointDownwash(Gaero,1),BoxPointDownwash(Gaero,2)); % G spline
            zout_div_dx_dw(Haero,1)=dTdx_Stickmodel(nodal_coordinates(H,1),nodal_coordinates(H,2),PHI(H,jshape),BoxPointDownwash(Haero,1),BoxPointDownwash(Haero,2)); % H spline           

            Qij=PHI(:,ishape)'*T_hFE_to_fp'*diag(ProjectedBoxArea)*inv(AIC)*(1i*RedFreq_lst(ii)/b*zout_dw+zout_div_dx_dw);   
            Q(ishape,jshape)=Qij;
        end
    end
    Q_lst(:,:,ii)=Q;
end

%% Write the intermediary output file:

save('IntermediaryOutput1.mat', 'Q_lst');