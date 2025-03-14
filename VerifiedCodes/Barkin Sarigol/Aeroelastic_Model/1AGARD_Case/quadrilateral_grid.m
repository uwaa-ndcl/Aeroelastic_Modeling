function [nodes] = quadrilateral_grid(x_corners, y_corners, n)
    % Generates an n x n grid inside a quadrilateral defined by corner points
    % x_corners, y_corners: Vectors of size 4 containing the (x,y) coordinates
    % n: Number of grid points along each dimension

    % Define a parametric grid (s, t) ranging from 0 to 1
    s = linspace(0, 1, n);
    t = linspace(0, 1, n);
    [S, T] = meshgrid(s, t);
    
    % Bilinear interpolation formula for a quadrilateral
    X = (1 - S) .* (1 - T) * x_corners(1) + S .* (1 - T) * x_corners(2) + ...
        S .* T * x_corners(3) + (1 - S) .* T * x_corners(4);
    
    Y = (1 - S) .* (1 - T) * y_corners(1) + S .* (1 - T) * y_corners(2) + ...
        S .* T * y_corners(3) + (1 - S) .* T * y_corners(4);
    
    X=reshape(X, [], 1);
    Y=reshape(Y, [], 1);

    nodes=[X,Y];


    % Plot the grid
    % figure; hold on; axis equal;
    % plot([x_corners x_corners(1)], [y_corners y_corners(1)], 'k-', 'LineWidth', 2); % Quadrilateral outline
    % plot(X, Y, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r'); % Grid points
    % title('Grid inside Quadrilateral');
    % xlabel('X'); ylabel('Y');
    % hold off;
end