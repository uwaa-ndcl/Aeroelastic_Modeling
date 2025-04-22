clear all

%% Construct System State Space (This Section Has Been Copied Over from the Large Wing Project and Adapted)

%% Inputs

filename = 'user_input.json';
fid = fopen(filename, 'r');
rawData = fread(fid, '*char')';  % Read the file as a string
fclose(fid);

% Decode the JSON string into a MATLAB structure
data = jsondecode(rawData);

ns = data.n_modes;
DoIUseGust = data.DoIUseGust;
Speeds=data.Speeds;
rho=data.rho;
StrainGauge_NodeNumbers=data.StrainGauge_NodeNumbers; %horizontal list of structural nodes corresponding the strain gauges
Accelerometer_NodeNumbers=data.Accelerometer_NodeNumbers;
StrainGauge_Orientations=data.StrainGauge_Orientations; %horizontal list with length equal to # of strain gauges
n_controlmodes=data.n_controlmodes;
lag_roots_lst=data.lag_roots_lst;   % list of beta_bar s
gust_lag_roots_lst=data.gust_lag_roots_lst; %list of gust_beta_bar s
RedFreq_lst=data.RedFreq_lst;  %%defined as wb/U
b=data.ref_semichord;


% Load intermediary outputs for the following variables needed: Pbarss, PbarGs, Pbarsc, Mss, Msc, Kss, Css
load('IntermediaryOutput2.mat');
load('MARGEI_Mass.mat')
%% Only for Marge I 


load('MARGEI_Mass.mat')
load('MARGEI_Stiffness.mat')  %<-change name for Marge II

b=b*0.0254;
%imperial_to_si_factor=4.44822*0.0254; % (lbf to Newtons)*(inch to m)
imperial_to_si_factor=(0.453592*0.0254)*(0.0254); % (lb*inch/s^2 to Newtons)*(inch to m)
Mss=Mss*imperial_to_si_factor;
Kss=Kss*imperial_to_si_factor;
Css=zeros(ns);
Msc=Msc*imperial_to_si_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Here, we construct Ap, Bc, BG for different speeds. 

nLagG=length(gust_lag_roots_lst);
nLag=length(lag_roots_lst);
nspeeds=length(Speeds);
nc=n_controlmodes; % get the number of columns of PHIc

% Prepare state space model matrices (eq. 1.38): [Ap], [Bc], {BG}
% for a list of speeds.


% Allocate storage space for arrays
% Total number of states
nstatesTotal=2*ns+ns*nLag+nLagG;
nstates=nstatesTotal;

% an array storing the [Ap] matrices for the different speeds
ApMatrices=zeros(nstates,nstates,nspeeds);
% an array storing the [Bc] matrices for the different speeds
BcMatrices=[];

if(nc > 0)
    BcMatrices=zeros(nstates,3*nc,nspeeds); %3*nc comes from [uc, s*uc, s^2*uc]^T
end

% an array storing the [BG] matrices for the different speeds
BGMatrices=[];

if(DoIUseGust == "Yes")
    BGMatrices=zeros(nstates,3,nspeeds);
end

% Array for storing values of highest frequencies (rad/sec) for which the
% math models are valid
OmegaMax=zeros(nspeeds,1);

for ispeed=1:nspeeds                                %     <------ Loop on Speeeds

    Uspeed=Speeds(ispeed);
        
    % Find the highest frequncy (rad / sec) for which the model is valid:
    OmegaMax(ispeed)=max(RedFreq_lst)*Uspeed/b;

    % create Mbarbarss,Cbarbarss, and Kbarbarss (equations 1.17, 1.18)
    % only the ns x ns part is needed
    Mbarbarss=zeros(ns,ns);
    Cbarbarss=zeros(ns,ns);
    Kbarbarss=zeros(ns,ns);
    %
    Mbarbarss(:,:)=Mss(:,:)-0.5*rho*b*b*Pbarss(:,:,3);
    Kbarbarss(:,:)=Kss(:,:)-0.5*rho*Uspeed*Uspeed*Pbarss(:,:,1);
    Cbarbarss(:,:)=Css(:,:)-0.5*rho*Uspeed*b*Pbarss(:,:,2);
    %
    Mbarbarsc=[];
    Cbarbarsc=[];
    Kbarbarsc=[];
    %
    if (nc > 0)      
        % create Mbarbarsc,Cbarbarsc, and Kbarbarsc (equations 1.17, 1.18)
        % only the NS x nc part is needed

        Mbarbarsc=zeros(ns,nc);
        Cbarbarsc=zeros(ns,nc);
        Kbarbarsc=zeros(ns,nc);
        %
        Mbarbarsc(:,:)=Msc(:,:)-0.5*rho*b*b*Pbarsc(:,:,3);
        Kbarbarsc(:,:)=-0.5*rho*Uspeed*Uspeed*Pbarsc(:,:,1);
        Cbarbarsc(:,:)=-0.5*rho*Uspeed*b*Pbarsc(:,:,2);
        % Note: no [Ksc] and [Csc] here.
        % end if loop on whether there are control surfaces (if nc ne 0)
    end
    
    % Find the inverse of Mbarbarss (required to convert eq. 1.28 to standard
    % state-space form)
    InvMbarbarss=inv(Mbarbarss);
    
    % Convert the Lag terms from their Roger fit values 
    % to their Laplace domain values (eq. 1.12)
    
    if nLag >0        
        beta=lag_roots_lst*Uspeed/b;   %List of betas    
    end

    betaG=[];
    if(DoIUseGust == "Yes")
      betaG=gust_lag_roots_lst*Uspeed/b;
    end

    % Create the aerodynamic lag equations for the structural and control
    % motions (eq. 1.23, 1.24)
    
    Ar=[];
    Brs=[];
    Brc=[];
    if nLag > 0   
        Ar=zeros(ns*nLag,ns*nLag);
        Brs=zeros(ns*nLag,ns);
        if nc > 0
            Brc=zeros(ns*nLag,nc);
        end
        
        for iLags=1:nLag
            istart=(iLags-1)*ns+1;
            iend=iLags*ns;
          
            for i=istart:iend
                Ar(i,i)=-beta(iLags);
            end
        
            for i=1:ns
                for j=1:ns
                    Brs(istart-1+i,j)=Pbarss(i,j,3+iLags);
                end    
            end
        
            if (nc > 0)
                for i=1:ns
                    for j=1:nc
                        Brc(istart-1+i,j)=Pbarsc(i,j,3+iLags);
                    end    
                end
            end
        
        end % end loop on iLags
    end % end nLag > 0 if statement
    
    
    

    ArG=[];
    BrG=[];
    if (DoIUseGust == "Yes")
    % Create the aerodynamic lag equations for the gust column (eqs. 1.25, 1.26)  
        if nLagG > 0
            ArG=zeros(nLagG,nLagG);
            BrG=zeros(nLagG,1);
            for i=1:nLagG
                ArG(i,i)=-betaG(i);
                BrG(i,1)=1;  
            end
        end
    end
    


    % Prepare the matrix [Ir} (eqs. 1.28 and 1.29)
    Ir=[];
    if nLag > 0    
        Ir=zeros(ns,ns*nLag);
        for iLags=1:nLag
            jstart=(iLags-1)*ns;
            for i=1:ns
                Ir(i,jstart+i)=1.;
            end  
        end
    end
    

    
    
    PGr=[];
    if DoIUseGust == "Yes"
    % Prepare the matrix PGr (eqs. 1.28,1.30)

        if nLagG >0 
            PGr=zeros(ns,nLagG);
            for i=1:ns
                for j=1:nLagG
                    PGr(i,j)=PbarGs(i,1,j+3);
                end
            end
        end
    end
    


    % Eqs. 1.32
    T21=-InvMbarbarss*Kbarbarss;
    T22=-InvMbarbarss*Cbarbarss;
    

    
    T2r=[];
    if nLag >0 
        T2r=0.5*rho*Uspeed*Uspeed*InvMbarbarss*Ir; 
    end



    T2Gr=[];
    if DoIUseGust == "Yes"
        if(nLag >0 )
            T2Gr=0.5*rho*Uspeed*InvMbarbarss*PGr;
        end
    end



    T2c=[];
    if (nc > 0)
        Temp=zeros(ns,3*nc);
        for i=1:ns
            for j=1:nc
                Temp(i,j)=Kbarbarsc(i,j);    
                Temp(i,j+nc)=Cbarbarsc(i,j);
                Temp(i,j+2*nc)=Mbarbarsc(i,j);
            end
        end
        T2c=-InvMbarbarss*Temp;
        clear Temp;
    end



    T2G=[];
    
    if DoIUseGust == "Yes"
        Temp=zeros(ns,3);
    
        Temp(:,1)=PbarGs(:,1,1); % PbarGs0
        Temp(:,2)=PbarGs(:,1,2); % PbarGs1
        Temp(:,3)=PbarGs(:,1,3); % PbarGs2
   
        T2G=zeros(ns,3);
        T2G=0.5*rho*Uspeed*InvMbarbarss*Temp;  %%%%%%% Multiply by b? ????????????????
        clear Temp;
    end
    
    
    %------------------------------------------------------------------------------
    % Prepare the plant state space matrices at a single speed 
    % to be stored in the arrays of plant matrices at all speeds
    % This code is derived from a code capable of generating state space
    % models for multiple velocities
    % 
    % Eqs. 1.34-1.38
    %
    % Total number of states = nstates
    %
    Ap=zeros(nstates,nstates);
    Bc=[];
    if (nc >0)
        Bc=zeros(nstates,3*nc);
    end
    
    BG=[];
    if(DoIUseGust == "Yes")
        BG=zeros(nstates,3);
    end
    

    % The Ap matrix ----------------------------------------------------
    
    
    for i=1:ns
        Ap(i,ns+i)=1;
    end
    
    for i=1:ns
        for j=1:ns
            Ap(i+ns,j)=T21(i,j);
            Ap(i+ns,j+ns)=T22(i,j);    
        end
    end
    

    if nLag > 0
        for i=1:ns
            for j=1:ns*nLag
                Ap(i+ns,2*ns+j)=T2r(i,j);    
            end
        end
    end


    if (DoIUseGust == "Yes") 
        for i=1:ns
            for j=1:nLagG
                Ap(i+ns,2*ns+nLag*ns+j)=T2Gr(i,j);    
            end
        end
    end
    

    if(nLag > 0) 
        for i=1:nLag*ns
            for j=1:ns
                Ap(2*ns+i,ns+j)=Brs(i,j);    
            end
        end
    
        for i=1:nLag*ns
            for j=1:nLag*ns
                Ap(2*ns+i,2*ns+j)=Ar(i,j);    
            end
        end
    end
    

    if (DoIUseGust == "Yes")
        if(nLagG > 0)
            for i=1:nLagG
                for j=1:nLagG
                    Ap(2*ns+nLag*ns+i,2*ns+nLag*ns+j)=ArG(i,j);   
                end
            end
        end
    end
    %
    % End of Ap matrix generation
    


    %%% [Bc] Generation if nc>0 -----------------------------------------
    if (nc > 0)
        Bc=zeros(nstates,3*nc);
        for i=1:ns
            for j=1:3*nc
                Bc(ns+i,j)=T2c(i,j);
            end
        end
    
        if(nLag >0)
            for i=1:nLag*ns
                for j=1:nc
                    Bc(2*ns+i,nc+j)=Brc(i,j);
                end
            end
        end

    end %End of [Bc] generation if statement


    
    %%% Generate [BG] if we are using gusts ------------------------------
    if DoIUseGust == "Yes"

        BG=zeros(nstates,3);
    
        for i=1:ns
            for j=1:3
                BG(ns+i,j)=T2G(i,j);    
             end
        end
    
        if(nLagG > 0)
            for i=1:nLagG
                BG(2*ns+nLag*ns+i,2)=BrG(i,1);  
            end
        end

    end %End of [BG] generation if statement
    



    % The ASE model - the state equations:--------------------------------
    %
    % {xdot]=[Ap]{x}+[Bc]{uc}+[BG]{uG}
    %
    % {uc} =   {qc}
    %          s{qc}
    %         s^2{qc}
    %
    % {uG} = wG
    %       s wG
    %      s^2 wG
    %
    % Store in the arrays that contain system matrices for all speeds
    % The Ap,Bc,BG matrices for the current speed ispeed are stored in the
    % arrays ApMatrices, BcMatrices, BGMatrices, where the third index identify
    % the speed case.
    %
    ApMatrices(:,:,ispeed)=Ap;
    %
    if(nc > 0)
        BcMatrices(:,:,ispeed)=Bc(:,:);    
    end
    
    if(DoIUseGust == "Yes")
        BGMatrices(:,:,ispeed)=BG(:,:);
    end

end % end loop over speeds

%%% The following commented part is unverified.
% %% Output Equation Generation        
% %%%We have strain gauges,accelerometers and hall effect sensor    
% 
% %Initialize the C_lst, Dcp_lst, DG_lst matrix lists for multiple speeds.
% num_strain=length(StrainGauge_NodeNumbers);
% num_accel=length(Accelerometer_NodeNumbers);
% dimC=1+num_accel+num_strain; %Hall effect+strain+accelerometer
% 
% C_lst=zeros(dimC,nstates,nspeeds);
% 
% Dcp_lst=[];
% if nc>0
%     dimBc=size(Bc,2); % number of columns of Bc
%     Dcp_lst=zeros(dimC,dimBc,speeds);
% end
% 
% DG_lst=[];
% if DoIUseGust == "Yes"
%     dimBG=size(BG,2); % number of columns of BG
%     DG_lst=zeros(dimC,dimBG,speeds);
% end
% 
% Tdisp=zeros(ns,nstates);
% Tvel=zeros(ns,nstates);
% 
% for ii=1:ns
%     Tdisp(ii,ii)=1.;
%     Tvel(ii,ii+ns)=1.;
% end
% 
% for ispeed=1:nspeeds
% 
%     %--------------------------------------------------------------------------------------------
%     %Hall effect sensor measures the structural mode with lowest frequency (rigid body rotation).
%     C_hall=zeros(1,nstates);
%     C_hall(1)=1;
% 
%     if nc>0 
%         Dcp_hall=zeros(1,dimBc);
%     end
% 
%     if DoIUseGust == "Yes"
%         dimBG=size(BG,2); % number of columns of Bc
%         DG_hall=zeros(1,dimBG);
%     end
% 
%     %-------------------------------------------------------------------------------------------
%     %Accelerometer outputs:      
%     % {Vector of accelerations}=[Cp]{xp}+[Dcp]{uc}+[DGp]{ug}
% 
%     Ap(:,:)=ApMatrices(:,:,ispeed);
% 
%     PHIsACC=zeros(num_accel,nstates);
%     PHIsACC(:,ns+1:2*ns)=PHIs(Accelerometer_NodeNumbers,:); %Structural mode shapes' contribution to nodes where the accelerometers are.
% 
%     CpACC=PHIsACC*Tvel*Ap;
% 
%     DcpACC=[];
%     if(nc >0 )
%         Bc(:,:)=BcMatrices(:,:,ispeed);    
%         DcpACC=PHIsACC*Tvel*Bc;
%     end
% 
%     DGpACC=[];
%     if(DoIUseGust == "Yes")
%         BG(:,:)=BGMatrices(:,:,ispeed);
%         DGpACC=PHIsACC*Tvel*BG;
%     end
% 
%     %-------------------------------------------------------------------------------------------
%     %Strain Gauge outputs:
% 
% 
% 
%     %-------------------------------------------------------------------------------------------
%     %Combining all matrices together for all sensors and for all speeds
% 
%     C_complete=[C_hall;CpACC;Cstrain];
%     C_lst(:,:,ispeed)=C_complete;
% 
%     if nc>0
%         Dc_complete=[Dcp_hall;DcpACC;Dcstrain];
%         Dcp_lst(:,:,ispeed)=Dc_complete;
%     end
% 
%     if(DoIUseGust == "Yes")
%         DG_complete=[DG_hall;DGpACC;DGstrain];
%         DG_lst(:,:,ispeed)=DG_complete;
%     end
% 
% end %End Loop over speeds



%% Save output file:

save('StateSpaceOutput.mat',  'ApMatrices'); %Return as well: 'C_lst',

if nc>0
    save('StateSpaceOutput.mat', 'BcMatrices' , '-append'); %Return as well: 'Dcp_lst',
end 

if(DoIUseGust == "Yes")
    save('StateSpaceOutput.mat', 'DG_lst','BGMatrices' , '-append');
end
