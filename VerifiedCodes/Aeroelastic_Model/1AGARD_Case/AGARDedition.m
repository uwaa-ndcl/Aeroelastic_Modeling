clc;
clear;
%% Read data

filename = 'user_input.json';
fid = fopen(filename, 'r');
rawData = fread(fid, '*char')';  % Read the file as a string
fclose(fid);

% Decode the JSON string into a MATLAB structure
data = jsondecode(rawData);

% Read all input variables from the string
Mach =data.Mach;
RedFreq_lst=data.RedFreq_lst;  %%defined as omega*b/U
b=data.ref_semispan;  
                                               
%%Add paths
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'DLM_helperfunctions'));
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'Spline_helperfunctions'));

%%%Get some data from the mesher:
[ProjectedBoxArea,BoxArea,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection] = DLM_Mesh_Updated();

%% This part emulates FEM.
n_modes=4;   % modes in order: wing twist, wing bend, tail roll, tail pitch

n_dim=30; % dimensional length of FEM nodes for each panel eg: lower wing as n_dim*n_dim structural points.

lower_wing_coordinates=quadrilateral_grid([2.75,0,2.25,3.7],[-1,-0.01,-0.01,-1],n_dim);
upper_wing_coordinates=quadrilateral_grid([2.75,0,2.25,3.7],[1,0.01,0.01,1],n_dim);

lower_tail_coordinates=quadrilateral_grid([4.3,3.1,4.4,4.65],[-1,-0.01,-0.01,-1],n_dim);
upper_tail_coordinates=quadrilateral_grid([4.3,3.1,4.4,4.65],[1,0.01,0.01,1],n_dim);
                                                                                   
nodal_coordinates= [lower_wing_coordinates;upper_wing_coordinates;lower_tail_coordinates;upper_tail_coordinates];
n_struct_nodes=length(nodal_coordinates);
n_aero_nodes=1200;

lowerwing = [1:n_dim^2];
upperwing =[n_dim^2+1:2*n_dim^2];
lowertail=[2*n_dim^2+1:3*n_dim^2];
uppertail=[3*n_dim^2+1:4*n_dim^2];


lowerwingaero=[1:300];
upperwingaero=[301:600];
lowertailaero=[601:900];
uppertailaero=[901:1200];

wing_coordinates=[lower_wing_coordinates;upper_wing_coordinates];
tail_coordinates=[lower_tail_coordinates;upper_tail_coordinates];

phi1=wing_coordinates(:,2)/b.*(wing_coordinates(:,1)/b-2.25.*abs(wing_coordinates(:,2)/b)-0.85);
phi2=wing_coordinates(:,2)/b.*abs(wing_coordinates(:,2)/b);
phi3=tail_coordinates(:,2)/b;

phi4=zeros(length(tail_coordinates),1);
for i=1:length(tail_coordinates)
    if tail_coordinates(i,2) == 0
        phi4(i,1)=(tail_coordinates(i,1)-3.75)/b;                                                                                           %%%%%<--------
    else
        phi4(i,1)=((tail_coordinates(i,1)-3.75)/b).*(tail_coordinates(i,2)/abs(tail_coordinates(i,2)));                                     %%%%%<--------
    end
end
PHI = [phi1,phi2,zeros(length(wing_coordinates),2);
       zeros(length(tail_coordinates),2),phi3,phi4];

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
%%%In this part, we add forces due to all modes (for each frequency)NC
Q_lst=zeros(n_modes,n_modes,length(RedFreq_lst));    %%list for different frequencies

for ii=1:length(RedFreq_lst)
    RedFreq=RedFreq_lst(ii);
    %Notice division by 'RedFreq/b' in the input of the function below instead of RedFreq.
    [AIC] = DLM_Find_AICs(Mach, RedFreq/b,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection);  
    Q=zeros(n_modes);

    for ishape=1:n_modes
        T_hFE_to_fp=zeros(n_aero_nodes,n_struct_nodes);
        T_hFE_to_fp(lowerwingaero,lowerwing)=T_Find(nodal_coordinates(lowerwing,1),nodal_coordinates(lowerwing,2),PHI(lowerwing,ishape),BoxSendPoint3(lowerwingaero,1),BoxSendPoint3(lowerwingaero,2)); % wing spline
        T_hFE_to_fp(upperwingaero,upperwing)=T_Find(nodal_coordinates(upperwing,1),nodal_coordinates(upperwing,2),PHI(upperwing,ishape),BoxSendPoint3(upperwingaero,1),BoxSendPoint3(upperwingaero,2)); % wing spline
        T_hFE_to_fp(lowertailaero,lowertail)=T_Find(nodal_coordinates(lowertail,1),nodal_coordinates(lowertail,2),PHI(lowertail,ishape),BoxSendPoint3(lowertailaero,1),BoxSendPoint3(lowertailaero,2)); % tail spline
        T_hFE_to_fp(uppertailaero,uppertail)=T_Find(nodal_coordinates(uppertail,1),nodal_coordinates(uppertail,2),PHI(uppertail,ishape),BoxSendPoint3(uppertailaero,1),BoxSendPoint3(uppertailaero,2)); % tail spline
            
        if RedFreq==0
            figure;
            hold on;
            scatter3(nodal_coordinates(:,1),nodal_coordinates(:,2),PHI(:,ishape),'.'); % Slightly transparent
            scatter3(BoxSendPoint3(:,1),BoxSendPoint3(:,2),T_hFE_to_fp*PHI(:,ishape),'.'); % Slightly transparent
            view(3);
            title('Spline at Force Points');
            hold off;
        end



        for jshape=1:n_modes
            T_hFE_to_dw=zeros(n_aero_nodes,n_struct_nodes);
            T_hFE_to_dw(lowerwingaero,lowerwing)=T_Find(nodal_coordinates(lowerwing,1),nodal_coordinates(lowerwing,2),PHI(lowerwing,jshape),BoxPointDownwash(lowerwingaero,1),BoxPointDownwash(lowerwingaero,2)); % wing spline
            T_hFE_to_dw(upperwingaero,upperwing)=T_Find(nodal_coordinates(upperwing,1),nodal_coordinates(upperwing,2),PHI(upperwing,jshape),BoxPointDownwash(upperwingaero,1),BoxPointDownwash(upperwingaero,2)); % wing spline
            T_hFE_to_dw(lowertailaero,lowertail)=T_Find(nodal_coordinates(lowertail,1),nodal_coordinates(lowertail,2),PHI(lowertail,jshape),BoxPointDownwash(lowertailaero,1),BoxPointDownwash(lowertailaero,2)); % tail spline
            T_hFE_to_dw(uppertailaero,uppertail)=T_Find(nodal_coordinates(uppertail,1),nodal_coordinates(uppertail,2),PHI(uppertail,jshape),BoxPointDownwash(uppertailaero,1),BoxPointDownwash(uppertailaero,2)); % tail spline
            zout_dw=T_hFE_to_dw*PHI(:,jshape);


            dT_hFE_to_dw_dx=zeros(n_aero_nodes,n_struct_nodes);
            dT_hFE_to_dw_dx(lowerwingaero,lowerwing)=dTdx_Find(nodal_coordinates(lowerwing,1),nodal_coordinates(lowerwing,2),PHI(lowerwing,jshape),BoxPointDownwash(lowerwingaero,1),BoxPointDownwash(lowerwingaero,2)); % wing spline
            dT_hFE_to_dw_dx(upperwingaero,upperwing)=dTdx_Find(nodal_coordinates(upperwing,1),nodal_coordinates(upperwing,2),PHI(upperwing,jshape),BoxPointDownwash(upperwingaero,1),BoxPointDownwash(upperwingaero,2)); % wing spline
            dT_hFE_to_dw_dx(lowertailaero,lowertail)=dTdx_Find(nodal_coordinates(lowertail,1),nodal_coordinates(lowertail,2),PHI(lowertail,jshape),BoxPointDownwash(lowertailaero,1),BoxPointDownwash(lowertailaero,2)); % tail spline
            dT_hFE_to_dw_dx(uppertailaero,uppertail)=dTdx_Find(nodal_coordinates(uppertail,1),nodal_coordinates(uppertail,2),PHI(uppertail,jshape),BoxPointDownwash(uppertailaero,1),BoxPointDownwash(uppertailaero,2)); % tail spline
            zout_div_dx_dw=dT_hFE_to_dw_dx*PHI(:,jshape);


            if ishape == jshape && RedFreq==0
                figure;
                hold on;
                scatter3(nodal_coordinates(:,1),nodal_coordinates(:,2),PHI(:,jshape),'.'); % Slightly transparent
                scatter3(BoxPointDownwash(:,1),BoxPointDownwash(:,2),zout_dw(:,1),'.'); % Slightly transparent
                view(3);
                title('Spline at Downwash Points');
                hold off;
                figure;
                scatter3(BoxPointDownwash(:,1),BoxPointDownwash(:,2),zout_div_dx_dw(:,1),'.'); % Slightly transparent
                view(3);
                title('Gradient Spline at Downwash Points');
            end 



            Qij=PHI(:,ishape)'*T_hFE_to_fp'*diag(ProjectedBoxArea)*inv(AIC)*(1i*RedFreq_lst(ii)/b*zout_dw+zout_div_dx_dw);   
            Q(ishape,jshape)=Qij;
        end
    end
    Q_lst(:,:,ii)=Q;
end

%% Normalize Q_lst with total surface area to match the results

%Q_lst=Q_lst/(sum(sum(diag(ProjectedBoxArea))));

Q_normal=real(Q_lst)+1i*imag(Q_lst)/1.5; 
Q_normal=-1*Q_normal/2; % Divide by 2 because you model both wing sides together, multiply by -1 based on definition

angles=angle(Q_normal)*180/pi;   %%% DIRECTLY COMPARE THIS TO TABLE THREE IN PCKFM PAPER
magnitudes=abs(Q_normal);    %%% DIRECTLY COMPARE THIS TO TABLE THREE IN PCKFM PAPER

% ishape=2;
% jshape=2;
% [AIC] = DLM_Find_AICs(Mach, 1.5,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection);
% 
% T_hFE_to_fp=zeros(n_aero_nodes,n_struct_nodes);
% T_hFE_to_fp(wingaero,wing)=T_Find(nodal_coordinates(wing,1),nodal_coordinates(wing,2),PHI(wing,ishape),BoxSendPoint3(wingaero,1),BoxSendPoint3(wingaero,2)); % wing spline
% T_hFE_to_fp(tailaero,tail)=T_Find(nodal_coordinates(tail,1),nodal_coordinates(tail,2),PHI(tail,ishape),BoxSendPoint3(tailaero,1),BoxSendPoint3(tailaero,2)); % tail spline
% 
% T_hFE_to_dw=zeros(n_aero_nodes,n_struct_nodes);
% T_hFE_to_dw(wingaero,wing)=T_Find(nodal_coordinates(wing,1),nodal_coordinates(wing,2),PHI(wing,ishape),BoxPointDownwash(wingaero,1),BoxPointDownwash(wingaero,2)); % wing spline
% T_hFE_to_dw(tailaero,tail)=T_Find(nodal_coordinates(tail,1),nodal_coordinates(tail,2),PHI(tail,ishape),BoxPointDownwash(tailaero,1),BoxPointDownwash(tailaero,2)); % tail spline
% zout_dw=T_hFE_to_dw*PHI(:,jshape);
% 
% dT_hFE_to_dw_dx=zeros(n_aero_nodes,n_struct_nodes);
% dT_hFE_to_dw_dx(wingaero,wing)=dTdx_Find(nodal_coordinates(wing,1),nodal_coordinates(wing,2),PHI(wing,ishape),BoxPointDownwash(wingaero,1),BoxPointDownwash(wingaero,2)); % wing spline
% dT_hFE_to_dw_dx(tailaero,tail)=dTdx_Find(nodal_coordinates(tail,1),nodal_coordinates(tail,2),PHI(tail,ishape),BoxPointDownwash(tailaero,1),BoxPointDownwash(tailaero,2)); % tail spline       
% zout_div_dx_dw=dT_hFE_to_dw_dx*PHI(:,jshape);
% 
% 
% Qij=PHI(:,ishape)'*T_hFE_to_fp'*diag(ProjectedBoxArea)*inv(AIC)*(1i*1.5/b*zout_dw+zout_div_dx_dw);   
% newvar=Qij/(sum(sum(diag(ProjectedBoxArea))));
