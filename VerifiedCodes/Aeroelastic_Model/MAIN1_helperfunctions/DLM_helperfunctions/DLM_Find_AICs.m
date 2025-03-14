function [AIC] = DLM_Find_AICs(Mach, RedFreq,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection)
    clear global;
    
    % Lifting Surface by the Doublet Lattice Method
    %
    %
    % Insert the meshing / geometry source code here before the DLM process
    %
    % ????????????????????????????
    %
    %
    % 2024/11/24
    %
    % The DLM Process
    %
    % At the end of the meshing / geometry process we have the coordinates of:
    %
    % Points 1&2 of each box: The ends of the vortex line at 0.25 chord
    % BoxSendPoint1(nBoxesTotal,3); 
    % BoxSendPoint2(nBoxesTotal,3);
    %
    % The force points (where the vorticies are)
    % BoxSendPoint3(nBoxesTotal,3);
    %
    % The downwash points, where normalwash BC are enforced
    % BoxPointDownwash(nBoxesTotal,3);
    %
    % BoxWidth=zeros(nBoxesTotal,1);
    % BoxPointChord=zeros(nBoxesTotal,1);
    %
    % Also:
    %
    % nPanels
    %
    % PanelPoints
    % PanelPoints=zeros(3,4,nPanels);
    %
    % nBoxesTotal - Total Number of boxes
    %
    % Number of chordwise divisions for each panel
    % NchordDivisions=zeros(nPanels,1);
    % Number of spanwise divisions for each panel
    % NspanDivisions=zeros(nPanels,1);
    % Number of boxes for each panel
    % NboxesPerPanel=zeros(nPanels,1);
    %
    % nBoxesMax=nBoxesTotal;
    %
    % Boxes(3,4,nBoxesMax);
    %
    % BoxArea (defined as a column vector)
    % BoxArea(nBoxesMax,1);
    %
    % ProjectedBoxArea(nBoxesMax,1);
    % Cosines of Dihedral angles of Boxes
    % CosDihedral=zeros(nBoxesMax,1);
    %
    % Unit normals to boxes 
    % (boxes that are vertical are consistent with the way normals are
    % defined for non verrtical boxes. That is, for verrtical boxes
    % the normal is in the -y direction)
    % normals(3,nBoxesMax)
    %
    % iBoxesTable=zeros(NSmax,NCmax,nPanels);
    % iBoxBelongsToPanel=zeros(nBoxesMax);
    % iBoxOnPanelChordwise=zeros(nBoxesMax);
    % iBoxOnPanelSpanwise=zeros(nBoxesMax);
    %
    im=1i; % The pure imaginary number
    %
    % The Doublet Lattice Process
    % ----------------------------
    
    
    %Run Geometry Meshing Code
    %[ProjectedBoxArea,BoxArea,BoxPointChord,nBoxesTotal,BoxPointDownwash,BoxSendPoint1,BoxSendPoint2,BoxSendPoint3,normals,FlowDirection] = DLM_Mesh_Updated();
    
    
    % ------------------------
    % Generate the AIC matrix
    % ------------------------
    %
    % Allocate memory:
    % 7 points per each couple of sending/receiving boxes. 
    % Each with three components
    %
    % On the sending box:
    Point1=zeros(3,1);
    Point2=zeros(3,1);
    Point3=zeros(3,1);
    % On the receiving box:
    Point4=zeros(3,1);
    Point5=zeros(3,1);
    Point6=zeros(3,1);
    Point7=zeros(3,1);
    %
    % Receiving boxes {upwash/U}=[AIC]{Gamma} Sending boxes ?????????????????
    % AICij represents the effect of the vortex or pressure at sending box j 
    % on receiving (upwash) box i
    %
    idBoxi=1;
    idBoxj=1;
    
    for idBoxi=1:nBoxesTotal % Loop over all boxes
        %
        % Loop over receiving boxes
        %
        % row idBoxi = downwash at downwash point on the idBoxi-th box
        % point 4 is the downwash (upwash) point
        %
        Point4(:)=BoxPointDownwash(idBoxi,:); % more compact command
        % for i=1:3
        %     Point4(i)=BoxPointDownwash(idBoxi,i);
        % end
        %
        % end points and center, load, point on vortex line of receiving box
        % Note that points 1,2,3 on the vortex line for all boxes are stored in
        % BoxSendPoint1,BoxSendPoint2,BoxSendPoint3. The geometry, then,
        % applies to both sending and receiving boxes.
        %
        % for i=1:3
        %     Point5(i)=BoxSendPoint1(idBoxi,i);
        %     Point6(i)=BoxSendPoint2(idBoxi,i);
        %     Point7(i)=BoxSendPoint3(idBoxi,i);
        % end
        Point5(:)=BoxSendPoint1(idBoxi,:);
        Point6(:)=BoxSendPoint2(idBoxi,:);
        Point7(:)=BoxSendPoint3(idBoxi,:);
    
    
        for idBoxj=1:nBoxesTotal % Loop over sending boxes
    
            %disp(['idBoxi: ', num2str(idBoxi)]);
            %disp(['idBoxj: ', num2str(idBoxj)]);
            %
            % Points 1,2 and 3 for the idBoxj-th box
            %
            % for ii=1:3
            %     Point1(ii)=BoxSendPoint1(idBoxj,ii);
            %     Point2(ii)=BoxSendPoint2(idBoxj,ii);
            %     Point3(ii)=BoxSendPoint3(idBoxj,ii);
            % end
            Point1(:)=BoxSendPoint1(idBoxj,:);
            Point2(:)=BoxSendPoint2(idBoxj,:);
            Point3(:)=BoxSendPoint3(idBoxj,:);
    
            %
            % Find induced velocity at point 4 (on the receiving box idBoxi)
            %
            % call the function that calculates AIC elements.
            % inputs to the function: 
            % Points 4,5,6,7 on the receiving boc
            % points 1,2,3 on the sending box
            % BoxPointChord(idBoxj) of the sending box
            % Mach number and reduced frequency:
            % from array data_aero given at input.
            % Note: Here we add a loop over the reduced frequencies in the
            % data_aero array. ???????????????????????
            % The AIC matrices for all required k values can be created element
            % by element or matrix by matrix. This should be added.
            % the code is written now for just one Mach number and one reduced
            % frequency given by:
            % data_aero=[ element1 element2 ];
            % Mach=data_aero(1)
            % RedFreq=data_aero(2)
            % ????????????????????????????????????????
            %
            %disp([Point1,Point2,Point3,Point4,Point5,Point6,Point7]);
            
            
            [Drs,indic1,indic2,D1rs,D2rs,A1,P1m,xsr,ybar,zbar,tanLambda,e,gamma_sr,delta1,delta2]=dlm_AICelement1972(Point1,Point2,Point3,Point4,Point5,...
                Point6,Point7,BoxPointChord(idBoxj),Mach,RedFreq);
            %Use the below if statement for verification, to print the value of
            %a critical variable at a particular idBoxi and idBoxj
            % if idBoxi == 188 && idBoxj == 2     %USED ONLY FOR VERIFICATION 
            %     disp(Drs);                   %USED ONLY FOR VERIFICATION
            %     disp(D1rs);                    %USED ONLY FOR VERIFICATION
            %     disp(D2rs);                   %USED ONLY FOR VERIFICATION
            %     disp(A1);                    %USED ONLY FOR VERIFICATION
            %     disp(P1m);                   %USED ONLY FOR VERIFICATION
            %     disp(xsr);                    %USED ONLY FOR VERIFICATION
            %     disp(ybar);                   %USED ONLY FOR VERIFICATION
            %     disp(zbar);                    %USED ONLY FOR VERIFICATION
            %     disp(gamma_sr);                %USED ONLY FOR VERIFICATION
            %     disp(tanLambda);                   %USED ONLY FOR VERIFICATION
            %     disp(e);                       %USED ONLY FOR VERIFICATION
            %     disp(Mach);                    %USED ONLY FOR VERIFICATION
            %     disp(RedFreq);                 %USED ONLY FOR VERIFICATION
            %     BoxPointChord(idBoxj);         %USED ONLY FOR VERIFICATION
            %     disp(delta1);                  %USED ONLY FOR VERIFICATION
            %     disp(delta2);                  %USED ONLY FOR VERIFICATION
            % 
            % end                               %USED ONLY FOR VERIFICATION
       
            AICelement(idBoxi,idBoxj)=Drs;
            % normal to receiving box idBoxi
            % 
            normalreceiving(1)=normals(1,idBoxi);
            normalreceiving(2)=normals(2,idBoxi);
            normalreceiving(3)=normals(3,idBoxi);
            %
            % The vortex lattice steady AIC connecting box idBoxi and idBoxj:
    
            AICelesteady(idBoxi,idBoxj)=0.5*BoxPointChord(idBoxj)*VLMInducedSpeed(Point1,...       %This is replacing the K20 and K10 we removed earlier in a more accurate way, ie using VLM.      
                Point2,Point4,FlowDirection,Mach,normalreceiving);                                 %This is the quasi-steady part of the solution (in phase with downwash). eg: A wing rotates                  
                                                                                                   %with velocity theta_dot about x=x0. At the point (xD,yD,zD), the normal velocity due to rotation is
                                                                                                   % =-theta_dot*(xD-x0). By definition of quasisteadiness, the corresponding load is the same as for a stationary surface so cambered that 
                                                                                                   % =-theta_dot*(xD-x0). This is why we can use VLM to model the quasisteady part.
            %Indic_array1(idBoxi,idBoxj)=indic1;    %USED ONLY FOR CODE VERIFICATION
            %Indic_array2(idBoxi,idBoxj)=indic2;    %USED ONLY FOR CODE VERIFICATION
            %D1rs_arr(idBoxi,idBoxj)=D1rs;          %USED ONLY FOR CODE VERIFICATION
        end % end of loop on idBoxj sending boxes
    end % end of loop idBoxi on receiving boxes
    %
    
    AIC=AICelesteady+AICelement;
    
    
    % %%%Now save the AIC into a .mat file. 
    % 
    % % Name the file using the Mach number and k (ie: omega/U).
    % filename = sprintf('M_%d_k_%d.mat', Mach, RedFreq);
    % 
    % % Save the array to the file
    % save(filename, 'AIC');
    
end

% %%%%%Below applies for a wing with two modes.
% %1) Plunge mode with amplitude =1.
% %2) Rotation with amplitude =1 radian about x=0.
% %%%%Mode 1 (subscript h)

% w_vec=1i*RedFreq*ones(nBoxesTotal,1);
% pressure=inv(AIC)*w_vec;
% LiftForces=diag(ProjectedBoxArea)*pressure; % We use projected area. This is the same as using the regular mesh area and multiplying it by cos(dihedral angle) to find the z-direction contribution of force (ie lift).
% CLh=sum(LiftForces)/sum(ProjectedBoxArea);
% 
% disp(abs(CLh));
% disp(angle(CLh)*180/pi);
% 
% %%%%Mode 2 (subscript a)
% 
% w_vec=ones(nBoxesTotal,1)+1i*RedFreq*BoxPointDownwash(:,1); 
% %This is a sum of the effect of freestream flow achieving a normal
% %component to the mesh due to rotation angle + downwash (ie a local pluge-like
% %effect) due to rotation speed.
% 
% pressure=inv(AIC)*w_vec;
% LiftMoment=diag(BoxSendPoint3(:,1))*diag(ProjectedBoxArea)*pressure; % moment about x=0. (BoxSendPoint3(:,1)) is the moment arm)
% CMa=sum(LiftMoment)/sum(ProjectedBoxArea);
% 
% LiftForces=diag(ProjectedBoxArea)*pressure;
% CLa=sum(LiftForces)/sum(ProjectedBoxArea);
% 
% disp(abs(CLa));
% disp(angle(CLa)*180/pi);
