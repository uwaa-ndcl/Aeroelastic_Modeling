%
function [T_matrix] = T_Find(xin,yin,zin,xout,yout) %%%One mode shape at a time
%
% Eli L. - January 20 2025
%
% Harder, R.L., and Demarais, R.N., "Interpolation Using Surface Splines", 
% Journal of Aircraft, February 1972, pp. 189-190
%
% Find surface spline coefficients over an area
% Make sure to have enough input points on the circumference of the area
% and even slightly beyond it.
%
% inputs: xi, yi corrdinates of npoints over the surface
% npoints
% xin(npoints) - the xi locations
% yin(xpoints) - the yi locations
%
% Shape of motion, normal to the surface, are given at
% the npoints input points column by column and stored in
% zin(npoints)
npoints=length(xin);

% Output:
% Gives T matrix where T*z_in=z_out
% A(xin)*coeffs=zin
% B(xin,xout)*coeffs=zout
%Hence B*inv(A)*zin=zout
%B*inv(A)=T


%
%
% A Measure of Small Distance (see eq. 12 in the Harder / Demarais paper)
%
EpsDistance=1.d-8;
%
%
% Generate the equations matrix A
%
% Allocate storage for the matrix [A]:
A=zeros(npoints+3,npoints+3);
%
% First three equations:
% Sum of Fi =0
% Sum of xin(i)*Fi=0.
% Sum of yin(i)*Fi=0.
%
for i=1:npoints % Equations 1 to 3
    A(1,i+3)=1.;
    A(2,i+3)=xin(i);
    A(3,i+3)=yin(i);
end
%
% The following npoints equations
%
% Motion at point j due to contributions from all points i
%
for j=1:npoints % equatyions 4 to 3+npoints
    A(j+3,1)=1.;
    A(j+3,2)=xin(j);
    A(j+3,3)=yin(j);
    for i=1:npoints
        % Rij squared
        RijSquared=(xin(i)-xin(j))^2+(yin(i)-yin(j))^2;
        if(RijSquared >= EpsDistance) ; % check for a very small distance
        A(j+3,3+i)=RijSquared*log(RijSquared);
        else
        % end if
        end
        % end loop on i
    end
    % end loop on j
end

%-------------------Do the same for Matrix B

% Allocate storage for the matrix [B]:
npoints_out=length(xout);

B=zeros(npoints_out+3,npoints+3);
%
% First three equations:
% Sum of Fi =0
% Sum of xin(i)*Fi=0.
% Sum of yin(i)*Fi=0.
%
for i=1:npoints % Equations 1 to 3
    B(1,i+3)=1.;
    B(2,i+3)=xin(i);
    B(3,i+3)=yin(i);
end
%
% The following npoints equations
%
% Motion at point j due to contributions from all points i
%
for j=1:npoints_out % equatyions 4 to 3+npoints
    B(j+3,1)=1.;
    B(j+3,2)=xout(j);
    B(j+3,3)=yout(j);
    for i=1:npoints
        % Rij squared
        RijSquared=(xin(i)-xout(j))^2+(yin(i)-yout(j))^2;
        if(RijSquared >= EpsDistance)  % check for a very small distance
        B(j+3,3+i)=RijSquared*log(RijSquared);
        else
        % end if
        end
        % end loop on i
    end
    % end loop on j
end

% Perform Singular Value Decomposition
[U, S, V] = svd(A);

% Invert S while handling small singular values
tol = 1e-10; % Threshold for numerical stability
S_inv = diag(1 ./ diag(S)); % Compute inverse of S
S_inv(diag(S) < tol) = 0; % Regularize small singular values

% Compute pseudo-inverse using SVD
A_pseudo_inv = V * S_inv * U';




T_matrix_unedited=B*A_pseudo_inv;
%%Because of the way in which T is defined, we get the matrix
%%in the shape (zout+3 by zin+3) such that
%%T_matrix_unedited*[0;0;0;zin]=[0;0;0;zout]
%%Very clearly, the first 3 equations (rows) are not needed since they don't
%%contribute to the relationship between zin and zout. Also, first three
%%columns of T_matrix_unedited multiply 0s so they are not needed too.
%%Hence:

T_matrix=T_matrix_unedited(4:length(xout)+3,4:length(zin)+3);


disp(cond(A));
end
%
% end function call
%

    
