function [new_im] = DrawLineBetweenTwoPoints(trajectory, im, color)
%DrawLineBetweenTwoPoints. Fast function to change the pixel value of a
%certain image to draw a line segment between first and last point.
if nargin < 4
    color = [0 0 0];
end
first = trajectory(1,1:2)
last = trajectory(2,1:2)
displacement_vector = last - first
distance = sqrt(displacement_vector * displacement_vector.')
number_of_waypoints = ceil(distance)
step_displacement = displacement_vector / number_of_waypoints;
new_im = im;
for i = 0:number_of_waypoints
    waypoint = first + step_displacement * i;
    new_im(floor(waypoint(2)), floor(waypoint(1)),:) = color;
    new_im(floor(waypoint(2)), ceil(waypoint(1)),:) = color;
    new_im(ceil(waypoint(2)), floor(waypoint(1)),:) = color;
    new_im(ceil(waypoint(2)), ceil(waypoint(1)),:) = color;
end
end

