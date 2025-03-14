function Vout = VLMInducedSpeed(Point1,Point2,Point3,FlowDirection,mach,...
             normalreceiving)
    
    % SW Feb 20 2012
    % EL 2-9-2019
    %
    % Warren F. Phillips: Mechanics of Flight 
    %
    % Finds the induced speed at {Point3} due to a horseshoe vortex
    % from {Point1} (left) to {Point2} (right)
    %
    %
    % Inputs:
    % Point1 contains the x,y,z coordinates of A
    % Point2 contains the x,y,z coordinates of B
    % Point3 contains the x,y,z coordinates of the point where the induced
    %        velocity is calculated
    % The attached vortex segment point1 - point2 (A - B)
    % FlowDirection - Direction of incoming flow in the x,y,z axes
    % FlowDirection=zeros(3,1); ?????
    % FlowDirection=[1. 0. 0.]'; ?????
    %
    % Output:
    % Vout is the normal speed at the receing downwash point
    %
    % V contains the components Vx,Vy,Vz of the induced speed at Point3
    %
    % Make sure that Point1,Point2, and Point3 are in column vector form
    
    
    
    % Step 1-----------------------------------------------------------------
    
    % The Glauert transformation of the x axis (to account for compressibility)
    % See Warren F. Phillips: Mechanics of Flight Equations 1.11.1 to 1.11.3.
    % Warren F. Phillips multiplies the the y coordinate by (1-mach^2) to get
    % to the incompressible Laplace equation. To the same effect, we divide the 
    % x coordinate by 1-mach^2 instead.
    apara=1/sqrt(1.0-mach^2); 
    Point1(1)=Point1(1)*apara;
    Point2(1)=Point2(1)*apara;
    Point3(1)=Point3(1)*apara;
    
    %Step 2------------------------------------------------------------------
    % Implement Equation 1.9.1 in Warren F. Phillips: Mechanics of Flight
    % Implement iff point 3 does not lie on the line connecting points 1 and 2
    % Otherwise raise warning.
    
    % vector from 1 to 3
    R1=Point3-Point1;
    % vector from 2 to 3
    R2=Point3-Point2;
    
    if ((abs(norm(R1,2)*norm(R2,2)+dot(R1,R2))) >= 0.00000001) % Check if point 3 is not on the line connecting points 1 and 2 
        V_bound=(norm(R1,2)+norm(R2,2))/norm(R1,2)/norm(R2,2)/(norm(R1,2)*norm(R2,2)+dot(R1,R2))*cross(R1,R2) ; %Part of Equation 1.9.1
        V_right_trailing= -cross(FlowDirection,R1)/norm(R1,2)/(norm(R1,2)-dot(FlowDirection,R1)) ; %Part of Equation 1.9.1
        V_left_trailing= cross(FlowDirection,R2)/norm(R2,2)/(norm(R2,2)-dot(FlowDirection,R2)) ; %Part of Equation 1.9.1
        V=(V_bound+V_right_trailing+V_left_trailing)/4./pi;   %Complete Equation 1.9.1
        Vout=dot(V,normalreceiving); %Part of downwash that is normal to receiving surface
    else
        % downwash point is on line 1-2
        error('Downwash point is on line 1-2. Not allowed.'); 
    end

end


% most modern (Desmarais) used.