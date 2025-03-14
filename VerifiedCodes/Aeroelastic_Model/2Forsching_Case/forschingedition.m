
%% Calculation of Energy from Forsching Graphs (Particularly, Fig. 3)

span=2*0.88; 
hingepoint=0.42; % in meters from leading edge
b=0.3; % semichord
angle=0.82*pi/180; % see Fig 3
Reduced_Freq=0.372; % see Fig 3

realcp=[0.4210912840572414	0.090136519;  % from Fig 3 ( the rows are in the form: [x location, real cp])
0.42632259251370025	 0.10469820491469961;
0.4384946069983153	0.076214442;
0.45407960728878227	0.052100463;
0.48725963865920485	0.027790709;
0.5214083770647354	 0.012840378768807776;
0.5464408125326776	0.005162758;
0.5840548207819368	-0.001674251;
0.6  0];    % due to Kutta condition

imcp=[0.41891188393428624, 0.008616030154327563; % see Fig 3 ( the rows are in the form: [x location, imaginary cp])
0.4373230922409501, 0.010384894388734795;
0.45453381694047373, 0.013336320318611768;
0.4868501529051986, 0.013919493634876613;
0.5464959817936135, 0.0096190882583031;
0.5865756347343716, 0.003550103122110815;
0.6 0]; % due to Kutta condition

real_distance=hingepoint-realcp(:,1);    %%sign here???
imag_distance=hingepoint-imcp(:,1);

real_displacement=real_distance*angle;
imag_displacement=imag_distance*angle;

real_energy_distribution=real_displacement.*realcp(:,2); 
imag_energy_distribution=imag_displacement.*imcp(:,2);

real_interp_extrap_func = @(xq) interp1(realcp(:,1), real_energy_distribution, xq, 'linear', 'extrap');
imag_interp_extrap_func = @(xq) interp1(imcp(:,1), imag_energy_distribution, xq, 'linear', 'extrap');

x_lower=0.42;
x_upper=0.6;

real_energy = integral(real_interp_extrap_func, x_lower, x_upper)*span;
imag_energy = integral(imag_interp_extrap_func, x_lower, x_upper)*span;

Modeenergy=real_energy+1i*imag_energy;

%energy = integral of Cp*displacement(x) dx *span 
% Here, we are assuming that pressure coefficient does not vary along span,
% a big assumption. 


%% Calculation of Energy via Spline+DLM:

nodal_coordinates1=[0.45,0.02,0; % These immidate the nodal coordinates coming from a FEM code. 
                    0.55,0.02,0; % nodal_coordinates1 is the for the upper half span of the wing
                    0.84,0.87,0;
                    0.94,0.87,0];

nodal_coordinates2=[0.45,-0.02,0; % nodal_coordinates2 is the for the lower half span of the wing (replacing the reflector plate with an actual half wing) 
                    0.55,-0.02,0;
                    0.84,-0.87,0;
                    0.94,-0.87,0];
phi1=zeros(4,1); % calculate the mode shape for the upper half of the wing (to emulate the mode shape that would be obtained from a FEM code)
for i=1:4
    A=2.1445;
    B=-1;
    C=-0.9007;
    x0=nodal_coordinates1(i,1);
    y0=nodal_coordinates1(i,2);
    dist_to_hinge=abs(A*x0+B*y0+C)/(sqrt(A^2+B^2));
    phi1(i,1)=dist_to_hinge*angle;
end

phi2=zeros(4,1); % calculate the mode shape for the lower half of the wing (to emulate the mode shape that would be obtained from a FEM code)
for i=1:4
    A=-2.1445;
    B=-1;
    C=0.9007;
    x0=nodal_coordinates2(i,1);
    y0=nodal_coordinates2(i,2);
    dist_to_hinge=abs(A*x0+B*y0+C)/(sqrt(A^2+B^2));
    phi2(i,1)=dist_to_hinge*angle;
end


upper_lst=[401:600,2001:2200]; % these are the DLM box indices (based on user input) that correspond to the upper control surface.
lower_lst=[601:800,2201:2400]; % these are the DLM box indices (based on user input) that correspond to the lower control surface.



%%Add paths
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'DLM_helperfunctions'));
addpath(fullfile(pwd, '..', 'MAIN1_helperfunctions', 'Spline_helperfunctions'));
%addpath('MAIN1_helperfunctions/DLM_helperfunctions');
%addpath('MAIN1_helperfunctions/Spline_helperfunctions');

%%%Get some data from the mesher:
[ProjectedBoxArea,BoxArea,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection] = DLM_Mesh_Updated();
[AIC] = DLM_Find_AICs(0, Reduced_Freq/b,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection);



zout_div_dx_dw=zeros(2400,1);

dT_hFE_to_dw_dx1=dTdx_Find(nodal_coordinates1(:,1),nodal_coordinates1(:,2),phi1(:,1),BoxPointDownwash(upper_lst,1),BoxPointDownwash(upper_lst,2));
zout_div_dx_dw(upper_lst)=dT_hFE_to_dw_dx1*phi1(:,1);   %<--- 0.013  but should be 0.0143??


dT_hFE_to_dw_dx2=dTdx_Find(nodal_coordinates2(:,1),nodal_coordinates2(:,2),phi2(:,1),BoxPointDownwash(lower_lst,1),BoxPointDownwash(lower_lst,2));
zout_div_dx_dw(lower_lst)=dT_hFE_to_dw_dx2*phi2(:,1);
 

zout_dw=zeros(2400,1);

T_hFE_to_dw1=T_Find(nodal_coordinates1(:,1),nodal_coordinates1(:,2),phi1(:,1),BoxPointDownwash(upper_lst,1),BoxPointDownwash(upper_lst,2));
zout_dw(upper_lst)=T_hFE_to_dw1*phi1(:,1);


T_hFE_to_dw2=T_Find(nodal_coordinates2(:,1),nodal_coordinates2(:,2),phi2(:,1),BoxPointDownwash(lower_lst,1),BoxPointDownwash(lower_lst,2));
zout_dw(lower_lst)=T_hFE_to_dw2*phi2(:,1);


T_hFE_to_fp=zeros(2400,8);
T_hFE_to_fp(upper_lst,1:4)=T_Find(nodal_coordinates1(:,1),nodal_coordinates1(:,2),phi1(:,1),BoxSendPoint3(upper_lst,1),BoxSendPoint3(upper_lst,2));
T_hFE_to_fp(lower_lst,5:8)=T_Find(nodal_coordinates2(:,1),nodal_coordinates2(:,2),phi2(:,1),BoxSendPoint3(lower_lst,1),BoxSendPoint3(lower_lst,2));


Qij=[phi1;phi2]'*T_hFE_to_fp'*diag(ProjectedBoxArea)*inv(AIC)*(1i*Reduced_Freq/b*zout_dw+zout_div_dx_dw);

disp('The mode energy (ie Q) calculated using Figure 3 pressure coefficient plot is');
disp(Modeenergy);
disp('whereas the Q calculated using the DLM+Spline is');
disp(Qij);


%% Recreate the pressure coefficient plot (Fig 3) using data from the DLM code for Further Verification

wgrad=zeros(2400,1);   % x direction-gradient of the mode shape (ie: w/U)  
wgrad(2001:2400)=-angle; %11,12
wgrad(401:800)=-angle; %3,4
% 5,6 front thin lines
% 7, 8 back thin lines

wdisp=zeros(841,1);  % the mode shape itself (ie the displacements of the control surfaces are nonzero, zero elsewhere)
for i=[401:600,2001:2200] %upper control surface: inner flap, outer flap, 433:504 back thin line                upper front thin line is 289:360
    A=2.1445;
    B=-1;
    C=-0.9007;
    x0=BoxPointDownwash(i,1);
    y0=BoxPointDownwash(i,2);
    dist_to_hinge=abs(A*x0+B*y0+C)/(sqrt(A^2+B^2));
    wdisp(i)=-1i*Reduced_Freq/b*angle*dist_to_hinge;
end
 
for i=[601:800,2201:2400]  %lower conrol surface: inner flap, outer flap, 505:576 back thin line            lower front thin line is 361 to 432
    A=-2.1445;
    B=-1;
    C=0.9007;
    x0=BoxPointDownwash(i,1);
    y0=BoxPointDownwash(i,2);
    dist_to_hinge=abs(A*x0+B*y0+C)/(sqrt(A^2+B^2));
    wdisp(i)=-1i*Reduced_Freq/b*angle*dist_to_hinge;
end


w=wgrad+wdisp;
Cp=inv(AIC)*w;

ymin=0.475; % to plot the pressure coefficients ONLY at the same spanwise location as Fig 3
ymax=0.5;

figure;
hold on; % plot the real part of the pressure coefficient, directly compare to Figure 3
scatter3(BoxPointDownwash(:,1),BoxPointDownwash(:,2),Cp); 
view(3);
ylim([ymin ymax])
hold off;

figure;
hold on; % plot the imaginary part of the pressure coefficient, directly compare to Figure 3
scatter3(BoxPointDownwash(:,1),BoxPointDownwash(:,2),imag(Cp));
view(3);
ylim([ymin ymax])
hold off;




