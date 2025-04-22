%% Note that the state space becomes unstable at the same speed and frequency where we observe flutter in the Ugw plot

clear all

%% U-g-omga solution Inputs
rho=1.225;

load('IntermediaryOutput1.mat');
n_modes=15;
Q=Q_lst(1:n_modes,1:n_modes,:); % Only look at the structural modes for flutter

load('MARGEI_Mass.mat');  % loads Mss and Msc
load('MARGEI_Stiffness.mat'); % load Kss and Ksc

%Mss=diag(diag(Mss));   %<---------------could be used to remove the
%nondiagonal components of the M and K matrix.
%Kss=diag(diag(Kss));


%%%Read the rest of the data from json
filename = 'user_input.json';
fid = fopen(filename, 'r');
rawData = fread(fid, '*char')';  % Read the file as a string
fclose(fid);

% Decode the JSON string into a MATLAB structure
data = jsondecode(rawData);

% Read all input variables from the string
RedFreq_lst=data.RedFreq_lst;  %%defined as wb/U
b=data.ref_semichord*0.0254;                      %Gust modeshape input is a column vector with dimensions: # of aero boxes on gust vanes only x 1

%imperial_to_si_factor=4.44822*0.0254; % (lbf to Newtons)*(inch to m)
imperial_to_si_factor=(0.453592*0.0254)*(0.0254); % (lb*inch/s^2 to Newtons)*(inch to m)

%If reduced frequency in the user input does not include 0, add it in.
if ~ismember(0, RedFreq_lst)
    RedFreq_lst = [RedFreq_lst; 0]; 
end

RedFreq_lst=sort(RedFreq_lst); 

U_lst=[];
g_lst=[];
omega_lst=[];
k_liste=[];


%%% Taking g=2*xi and xi= 0.01;
%struct_damping=diag(diag(Kss))*imperial_to_si_factor*0.01*2*1i;   %This does not make much difference to frequencies.
struct_damping=zeros(n_modes);

%%Play around with Kss:
%Kss=diag(diag(Kss));
%Kss(5,5)=Kss(5,5)/2.3; %% First Fuselage bending
%Mss=diag(diag(Mss));

%% Uninterpolated Version
% for i=1:length(RedFreq_lst)   %% Start from 2 to avoid k=0?
% 
%     R=-Mss*imperial_to_si_factor -0.5*rho*b^2/(RedFreq_lst(i))^2*Q(:,:,i);
%     S=struct_damping+Kss*imperial_to_si_factor;
%     lambdas=eig(-R, S);
%     disp(lambdas)
%     for j=1:length(lambdas)
%         lambda=lambdas(j);
%         if real(lambda)>0
%             omega=sqrt(1/real(lambda));
%             g=imag(lambda)/real(lambda);
%             U=omega*b/RedFreq_lst(i);
%             U_lst=[U_lst,U];
%             g_lst=[g_lst,g];
%             omega_lst=[omega_lst,omega];
%             k_liste=[k_liste,RedFreq_lst(i)];
%         end
%     end
% end

%% Interpolated version
Fine_RedFreq_lst=linspace(min(RedFreq_lst),max(RedFreq_lst),2000);
Fine_Q=zeros(n_modes,n_modes,length(Fine_RedFreq_lst));

for i=1:n_modes
    for j=1:n_modes
        Q_row_col=reshape(Q(i,j,:),[length(RedFreq_lst),1]);
        Fine_Q(i,j,:) = interp1(RedFreq_lst, Q_row_col, Fine_RedFreq_lst);
    end
end

for i=1:length(Fine_RedFreq_lst) %% Start from 2 to avoid k=0?

    R=-Mss*imperial_to_si_factor-0.5*rho*b^2/(Fine_RedFreq_lst(i))^2*Fine_Q(:,:,i);
    S=struct_damping+Kss*imperial_to_si_factor;
    lambdas=eig(-R, S);
    disp(lambdas)
    for j=1:length(lambdas)
        lambda=lambdas(j);
        if real(lambda)>0
            omega=sqrt(1/real(lambda));
            g=imag(lambda)/real(lambda);
            U=omega*b/Fine_RedFreq_lst(i);
            U_lst=[U_lst,U];
            g_lst=[g_lst,g];
            omega_lst=[omega_lst,omega];
            k_liste=[k_liste,Fine_RedFreq_lst(i)];
        end
    end
end




%% Plotting
figure;
scatter(U_lst,g_lst,'Marker', '.');
xlabel('U (m/s)');
ylabel('g ');
xlim([0 100]);
ylim([-3 0.5]);

figure;
scatter(U_lst,omega_lst,'Marker', '.');
xlabel('U (m/s)');
ylabel('omega (rad/s)');
xlim([0 100]);
ylim([0 200]);
            

%Eigenvectors not exactly [1 0 0 0 .... 0]^T or [0 1 0 0 0 ....]^T because
%shapes of the structural modes are different to those of aerodynamics but
% you can see heavy contributions to eigenvectors (aeroelastic modes) from some pure structural
% modes as opposed to others. eig(-R, S)

% Normally you would not read the frequencies of modes from a Ugw graph
% because its'solutions are not physical but since damping is slow at
% around 0.02, the results are basically SHM which is only the case when
% Ugw g and w mean anything.



% %%% FRF %%% to see natural frequencies:
% q_lst=[];
% %%%%1) Pick a U value: eg 20 m/s
% Ums=20;
% 
% %%%%2) loop over k values and solve for q , (R+onedivomega^2*S)q=0.
% for i=1:length(RedFreq_lst)
% 
%     R=-Mss-0.5*rho*b^2/RedFreq_lst(i)^2*Q(:,:,i);
%     S=struct_damping+Kss;
% 
%     onedivomega= b/RedFreq_lst(i)*Ums;
% 
%     q=null(R+onedivomega^2+S);
%     q_lst=[q_lst,q];
% end
% 
% %%%3) Now plot some component of q against omega (defined as k*Ums/b)
% 
% figure;
% scatter(q_lst(1,:),RedFreq_lst*Ums/b);