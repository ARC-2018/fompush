function [drawn_image] = DrawTrajectoryOnImage(trajectory, im, color)
%DrawTrajectoryOnImage Gets the x and y values of the trajectory waypoints
% and draws it on image. It returns the modified image. If not set, color
% will be automatically set to blue
if nargin < 4
    color = [0 0 255];
end
% x_waypoints = trajectory(:,2);
% y_waypoints = trajectory(:,1);
drawn_image = im;     % creates a merged image
N = size(trajectory, 1);
for i = 2:N
    line_trajectory = trajectory((i-1):i,:);
    drawn_image = DrawLineBetweenTwoPoints(line_trajectory, drawn_image, color);
%     drawn_image(floor(x_waypoints(i)), floor(y_waypoints(i)),:) = color;
%     drawn_image(floor(x_waypoints(i)), ceil(y_waypoints(i)),:) = color;
%     drawn_image(ceil(x_waypoints(i)), floor(y_waypoints(i)),:) = color;
%     drawn_image(ceil(x_waypoints(i)), ceil(y_waypoints(i)),:) = color;
end
end

