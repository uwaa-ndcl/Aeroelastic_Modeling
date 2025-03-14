function [Drs,indicator1,indicator2,D1rs,D2rs,A1,P1m,xsr,ybar,zbar,tanLambda,e,gamma_sr,delta1,delta2] = dlm_AICelement1972(Point1,Point2,Point3,Point4,Point5,...
    Point6,Point7,BoxPointChord,Ma,k) %k = omega/U 
    %Note: This function normally only needs to return Drs. The rest of the
    %returned variables are used for verification in DLM_Find_AICs.m file


    % Inputs:
    % Point1 contains the x,y,z coordinates of left terminus of sending box
    % Point2 contains the x,y,z coordinates of right terminus of sending box
    % Point3 contains the x,y,z coordinates of middle terminus of sending box
    % Point4 contains the x,y,z coordinates of the point where the induced
    %        velocity is calculated (receiving box)
    % Point5 contains the x,y,z coordinates of left terminus of receiving box
    % Point6 contains the x,y,z coordinates of right terminus of receiving box
    % Point7 contains the x,y,z coordinates of middle terminus of receiving box
    % Point5,6,7 are used to calculate the dihedral angle of the receiving box 

    e = 0.5*sqrt((Point2(3)-Point1(3))^2+(Point2(2)-Point1(2))^2);
    
    % cartesian coordinates of receiving points relative to sending points
    %%%% Check this, not sure about point 3!!!!!!
    xsr = Point4(1)-Point3(1);   
    ysr = Point4(2)-Point3(2);  
    zsr = Point4(3)-Point3(3);
    %disp('start');
    %disp(xsr);
    %disp(Point3);
    %disp(Point4);
    %disp('end');
    % dihedral angle gamma = arctan(dz/dy) and sweep angle lambda = arctan(dx/dy)
    sinGamma = (Point2(3)-Point1(3)) / (2.0 * e);  %!!!!!!!
    cosGamma = (Point2(2)-Point1(2)) / (2.0 * e);  %!!!!!!!
    tanLambda = (Point2(1)-Point1(1)) / (2.0 * e); %!!!!!!!
    
    GAMAR=atan((Point7(3)-Point5(3))/(Point7(2)-Point5(2)));  %!!!!!!! CHECK THESE LINES.
    GAMAS=atan((Point3(3)-Point1(3))/(Point3(2)-Point1(2)));  %!!!!!!!   
    gamma_sr=GAMAS-GAMAR;
        

    % local coordinates of receiving point relative to sending point
    ybar = ysr * cosGamma + zsr * sinGamma;
    zbar = zsr * cosGamma - ysr * sinGamma;

    % pre-calculate some values which will be used a couple of times
    ybar2 = ybar^2.0;
    zbar2 = zbar^2.0;
    e2 = e^2.0;

    %Calculate Deltas in Rodden 1972 equation 30
    if ybar2+zbar2-e2>0
        delta1=1;
        delta2=0;
    elseif ybar2+zbar2-e2 == 0
        delta1=0;
        delta2=0.5;
    else
        delta1=1;
        delta2=1;
    end

    % Calculate Fparabolic
    ratio = 2.0 * e * abs(zbar) / (ybar2 + zbar2 - e2); 
   
    if abs(zbar)/e <= 0.001    %Condition 1 planar (ie i0)
        Fparabolic = 2.0*e/(ybar2 - e2);
        indicator1=0;              %Indicator1 and 2 variables are used purely for verification
    elseif abs(ratio) <= 0.3 && abs(zbar)/e > 0.001 % Condition 2, co-planar / close-by (ie ia)
      
        funny_series=0;
        for n = 2:7
            funny_series=funny_series + (-1.0)^n / (2.0 * n - 1.0) * ratio^(2.0 * n - 4.0);
            % Rodden 1971, eq 33, Rodden 1972, eq 31b and Rodden 1998, eq 25
        end
        epsilon = 4.0 * e^4.0 / (ybar2 + zbar2 - e2)^2.0 * funny_series;
        Fparabolic = delta1 * 2.0 * e / (ybar2 + zbar2 - e2) * (1.0 - epsilon * zbar2 / e2) +delta2*pi/abs(zbar);
        indicator1=1;
    elseif abs(ratio) > 0.3 && abs(zbar)/e > 0.001 % Condition 3, rest / further away (ie ir)
        
        epsilon=e2/zbar2*(1-((ybar2+zbar2-e2)/(2*e*abs(zbar)))*atan((2*e*abs(zbar))/(ybar2+zbar2-e2)));
        Fparabolic = delta1 * 2.0 * e / (ybar2 + zbar2 - e2) * (1.0 - epsilon * zbar2 / e2) +delta2*pi/abs(zbar);
        indicator1=2;
    else
        error('Something went wrong, check dlmAICelement1972 code');
    end

    % call the kernel function
    [P1m, P2m] = Kernel_Corrected(Ma,k,xsr,ybar,zbar, gamma_sr, tanLambda, -e);
    [P1p, P2p] = Kernel_Corrected(Ma,k,xsr,ybar,zbar, gamma_sr, tanLambda, e);
    [P1s, P2s] = Kernel_Corrected(Ma,k,xsr,ybar,zbar, gamma_sr, tanLambda, 0);
    %disp(P1s);
    % define terms used in the parabolic approximation, Nastran incro.f
    A1 = (P1m - 2.0 * P1s + P1p) / (2.0 * e2);  % Rodden 1971, eq 28
    B1 = (P1p - P1m) / (2.0 * e);  % Rodden 1971, eq 29
    C1 = P1s;  % Rodden 1971, eq 30

    A2 = (P2m - 2.0 * P2s + P2p) / (2.0 * e2);  % Rodden 1971, eq 37
    B2 = (P2p - P2m) / (2.0 * e);  % Rodden 1971, eq 38
    C2 = P2s;  % Rodden 1971, eq 39


    % The "planar" part, D1rs
    % -----------------
    %  normalwash matrix, Rodden 1971, eq 34
    D1rs = BoxPointChord / (pi * 8.0) *(((ybar2 - zbar2) * A1 + ybar * B1 + C1) * Fparabolic ...
        + (0.5 * B1 + ybar * A1) * log(((ybar - e)^2.0 + zbar2) / ((ybar + e)^2.0 + zbar2)) ...
        + 2.0 * e * A1);
    
    % The "nonplanar" part, D2rs
    % --------------------
    
    %Distinguish between two conditions:
    if abs(zbar)/e <= 0.001
        D2rs=0;   %%%%This is because the Kernel function gives P2=0 when zbar=0 is inputted. Just try it. 
                  %%%%D2rs is the integral of (A2*eta2+B2*eta+C2)/r^2 see
                  %%%%equation 33 in the 1972 paper.
                  %%%%Looking at lines 88 to 90, notice A2, B2, C2 =0 when
                  %%%%P2=0. Hence D2rs becomes 0.
        indicator2=0;
    else
        if abs(1/ratio) <= 0.1
            D2rs = BoxPointChord / (16.0 * pi * zbar2) * (((ybar2 + zbar2) * A2 + ybar * B2 + C2) * Fparabolic ...
                + 1.0 / ((ybar + e)^2.0 + zbar2) * (((ybar2 + zbar2) * ybar + (ybar2 - zbar2) * e)* A2 ...
                + (ybar2 + zbar2 + ybar * e) * B2 + (ybar + e) * C2) ...
                - 1.0 / ((ybar - e)^2.0 + zbar2)* (((ybar2 + zbar2) * ybar - (ybar2 - zbar2) * e)* A2 ...
                + (ybar2 + zbar2 - ybar * e) * B2 + (ybar - e) * C2)); 
            indicator2=1;
        else
            if ratio < 0.3
                epsilon=0;
                funny_series=0;
            
                for n = 2:7
                    funny_series=funny_series + (-1.0)^n / (2.0 * n - 1.0) * ratio^(2.0 * n - 4.0);
                    % Rodden 1971, eq 33, Rodden 1972, eq 31b and Rodden 1998, eq 25
                end
                epsilon = 4.0 * e^4.0 / (ybar2 + zbar2 - e2)^2.0 * funny_series;
                indicator2=2;
            else
                epsilon=e2/zbar2*(1-((ybar2+zbar2-e2)/(2*e*abs(zbar)))*atan((2*e*abs(zbar))/(ybar2+zbar2-e2)));
                indicator2=3;
            end
         
            DELTA=e2/zbar2*(1-delta1-delta2*pi/abs(zbar)*((ybar2+zbar2-e2)/(2*e)));
    
            D2rs = BoxPointChord * e / (8.0 * pi * (ybar2 + zbar2 - e2)) * ((2.0 * (ybar2 ...
                + zbar2 + e2) * (e2 * A2 + C2) + 4.0 * ybar * e2 * B2)/ (((ybar + e)^2.0 + zbar2) * ((ybar - e)^2.0 + zbar2)) ...
                - (delta1*epsilon+DELTA) / e2 * ((ybar2 + zbar2) * A2 + ybar * B2 + C2));
        end
    end
    % add planar and non-planar parts, % Rodden eq 22
    % the steady part D0 has already been subtracted inside the kernel function
    Drs = D1rs + D2rs;
    %disp([D1rs,D2rs]);
end





















   % 
    % elseif (abs(1.0 / ratio) > 0.1) && (abs(zbar) / e > 0.001)  % ie condition ic, Rodden 19671 eq 41   %%%CHECK THIS!!!
    %     if ratio <= 0.3
    %         epsilon=0;
    %         funny_series=0;
    % 
    %         for n = 2:7
    %             funny_series=funny_series + (-1.0)^n / (2.0 * n - 1.0) * ratio^(2.0 * n - 4.0);
    %             % Rodden 1971, eq 33, Rodden 1972, eq 31b and Rodden 1998, eq 25
    %         end
    %         epsilon = 4.0 * e^4.0 / (ybar2 + zbar2 - e2)^2.0 * funny_series;
    % 
    %     elseif ratio > 0.3
    %         Fparabolic = 1.0/abs(zbar)*arctan2((2.0*e*abs(zbar)), (ybar2+zbar2-e2));
    %         epsilon=(1.0 - Fparabolic * (ybar2 + zbar2 - e2) / (2.0 * e)) / zbar2 * e2;
    % 
    %     else
    %         error('Something went wrong, check dlmAICelement1972 code');
    %     end
    %     D2rs = chord * e / (8.0 * pi * (ybar2 + zbar2 - e2)) * ((2.0 * (ybar2 ...
    %         + zbar2 + e2) * (e2 * A2 + C2) + 4.0 * ybar * e2 * B2)/ (((ybar + e)^2.0 + zbar2) * ((ybar - e)^2.0 + zbar2)) ...
    %         - epsilon / e2 * ((ybar2 + zbar2) * A2 + ybar * B2 + C2));
    % else
    %     error('Something went wrong, check dlmAICelement1972 code');
    % end