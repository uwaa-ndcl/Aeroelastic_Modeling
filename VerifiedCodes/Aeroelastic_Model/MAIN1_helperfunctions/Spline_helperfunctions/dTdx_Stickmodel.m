
%Verification of this code has been done by comparing the hand calculated
%gradients on Marge I, panel BC for mode 9 to those provided by this code. 

%Example usage: zout_div_dx_dw(BCaero,1)=dTdx_Stickmodel(nodal_coordinates(BC,1),nodal_coordinates(BC,2),PHI(BC,jshape),BoxPointDownwash(BCaero,1),BoxPointDownwash(BCaero,2));


function dzdx_out = dTdx_Stickmodel(x_in, y_in, z_in, x_out, y_out)
    % Create scattered interpolant with extrapolation
    F = scatteredInterpolant(x_in, y_in, z_in, 'linear', 'linear');
    
    % Small tolerance for grouping points at the same y-location
    % Because we have a stick model, dw/dx = constant at a given y location
    
    tol = 1e-8;
    
    % Initialize dzdx values as NaN (ignored if no valid gradient is found)
    dzdx_vals = NaN(size(z_in));
    
    % Loop through unique y values (considering tolerance)
    unique_y = unique(y_in);
    for i = 1:length(unique_y)
        y_val = unique_y(i);
        idx = abs(y_in - y_val) < tol;
        x_subset = x_in(idx);
        z_subset = z_in(idx);
        
        % calculate dw/dx from points where there is more than 1 point at
        % the same y (spanwise) location.
        if length(x_subset) > 1
            % Compute dz/dx using finite differences when multiple x points exist
            [sorted_x, sort_idx] = sort(x_subset);
            %disp(y_val);
            sorted_z = z_subset(sort_idx);
            format long
            %disp(x_subset)
            format long
            %disp(z_subset)
            dzdx_subset = gradient(sorted_z, sorted_x); 
            %disp(dzdx_subset);
            dzdx_vals(idx) = dzdx_subset;
            % The Matlab gradient fcn works such that:
            % If there are three points at a spanwise location the gradient
            % between points 1 and 2 is allocated to be the dzdx value at
            % point 1 and the gradient between points 2 and 3 is allocated
            % to be the gradient at point 3. Gradient at point 2 is halway
            % between gradient at point 1 and 3.
        end
    end
    %disp(dzdx_vals)
    % Create scattered interpolant for dz/dx (ignoring NaNs)
    valid_idx = ~isnan(dzdx_vals);
    F_x = scatteredInterpolant(x_in(valid_idx), y_in(valid_idx), dzdx_vals(valid_idx), 'linear', 'linear');
    %disp(x_in(valid_idx))
    %disp(x_in)
    % Interpolate dz/dx at (x_out, y_out)
    dzdx_out = F_x(x_out, y_out);
end


