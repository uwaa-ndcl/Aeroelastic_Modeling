function [P1,P2] = Kernel_Corrected(M,k,xbar,ybar,zbar, gamma_sr, tanLambda, ebar)
    
    % % This is the function that calculates "the" kernel function(s) of the DLM.
    % % K1,2 are reformulated in Rodden 1972 compared to Rodden 1969 and include new
    % % conditions, e.g. for co-planar panels.
    % % Note that the signs of K1,2 are switched in Rodden 1971, in this implementation
    % % we stay with the 1969 convention 
    % % from the VLM. Also, we directly subtract the steady parts K10 and K20, as the
    % % steady contribution will be added later from the VLM.

    r1 = ((ybar - ebar)^2.0 + zbar^2.0)^0.5;  % Rodden 1971, eq 4
    beta2 = 1.0 - (M^2.0);  % Rodden 1971, eq 9
    R = ((xbar - ebar * tanLambda)^2.0 + beta2 * r1^2.0)^0.5 ; % Rodden 1971, eq 10
    u1 = (M * R - xbar + ebar * tanLambda) / (beta2 * r1) ;  % Rodden 1971, eq 11
    k1 = k * r1 ; % Rodden 1971, eq 12 with k = w/U
    j = 1i ;  % imaginary number
    ejku = exp(-j * k1 * u1)  ;% pre-multiplication
    
    % direction cosine matrices
    T1 = cos(gamma_sr) ; % Rodden 1971, eq 5
    T2 = zbar * (zbar * cos(gamma_sr) + (ybar - ebar) * sin(gamma_sr)) ; % Rodden 1971, eq 21a: T2_new = T2_old*r1^2

    % % Approximation of intergrals I1,2, Rodden 1971, eq 13+14
    % [I1, I2] = I1I2TG(u1,k1);                                                           


    %Evaluate K1 and K2 considering the possible singularity when r1=0.
    if r1 == 0 && xbar >= 0.0
        K1= -2.0 ;                   %You get this by using the expression K1= I1 -ejku*M*r1.... The -ejku*M*r1 part equals 0 when r1=0. 
        K2= +4.0 ;                   %But the I1 evaluated using the I1I2TG.m function gives +2 when k1==0 and u1<0 (ie r1==0 and u1<0 since k1=w*r1/U) and gives 0 when r1==0 and u1>0. 
    elseif r1==0 && xbar< 0.0        %Using the line 14 expression for u1 and substituing in the line 13 expression for R, plotting u1=f(r1) and playing around with the value of x0, you see that
        K1 = 0.0 ;                   % approaching r1=0 from the positive r1 (because r1>=0 always by definition) yields u1= negative infinity at r1=0 when x0>0. So, I1I2TG(u1=-100000000000,k1=0) =2.
        K2 = 0.0 ; 
    else
        % Approximation of intergrals I1,2, Rodden 1971, eq 13+14
        [I1, I2] = I1I2TG(u1,k1); 
        % Formulation of K1,2 by Landahl, Rodden 1971, eq 7+8
        K1 = -I1 - ejku * M * r1 / R / (1 + u1^2.0)^0.5 ;
        K2 = 3.0 * I2 + j * k1 * ejku * (M^2.0) * (r1^2.0) / (R^2.0) / (1.0 + u1^2.0)^0.5 + ejku * M * r1 * ((1.0 + u1^2.0) * beta2 * r1^2.0 / R^2.0 + 2.0 + M * r1 * u1 / R)/(R*(1.0 + u1^2.0)^1.5) ;
    end

    % This is the analytical solution for K1,2 at k=0.0, Rodden 1971, eq 15+16
    K10 = -1.0 - (xbar - ebar * tanLambda) / R ;
    K20 = 2.0 + (xbar - ebar * tanLambda) * (2.0 + beta2 * r1^2.0 / R^2.0) / R ; 

    % Rodden 1971, eq 27b, check: -K1*exp(-j*k*xbar)*T1
    P1 = -(K1 * exp(-j * k * (xbar - ebar * tanLambda)) - K10) * T1 ;
    % Rodden 1971, eq 36b, check: -K2*exp(-j*k*xbar)*T2/r1**2.0
    P2 = -(K2 * exp(-j * k * (xbar - ebar * tanLambda)) - K20) * T2 ;

end
