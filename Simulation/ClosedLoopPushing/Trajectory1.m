function [x_star] = Trajectory1(t)

import Models.QSPusherSlider
u_thrust = 0.05;
c = QSPusherSlider.c;
px = QSPusherSlider.a / 2.0;
x_star = zeros(4, length(t));

t_1 = 2;
t_2 = 4;
t_3 = 10;
circle = @(d, A, B, C, t, t_0, x_0)(x_0 -[(sin(-(d*A)*t) * B + cos(-(d*A)*t) * C - sin(-(d*A)*t_0) * B - cos(-(d*A)*t_0) * C) / d;...
                                         (-cos(-(d*A)*t) * B + sin(-(d*A)*t) * C + cos(-(d*A)*t_0) * B - sin(-(d*A)*t_0) * C) / d;...
                                          (d*A) * (t - t_0);
                                          0]);
line = @(t, t_0, x_0)(x_0 + [u_thrust * (t - t_0) * cos(x_0(3)); u_thrust * (t - t_0) * sin(x_0(3)); 0; 0]);
d1 = 0.0249;
A1 = u_thrust / (c^2 + px^2 + d1^2);
B1 = (c^2 + px^2);
C1 = px*d1;
x0 = [0;0;0;d1];
x_1 = circle(d1, A1, B1, C1, t_1, 0, x0);
x_1(4) = 0;
d2 = -0.02;
A2 = u_thrust / (c^2 + px^2 + d2^2);
B2 = (c^2 + px^2);
C2 = px*d2;
x_2 = line(t_2, t_1, x_1);
x_2(4) = d2;
d3 = 0.023;
A3 = u_thrust / (c^2 + px^2 + d3^2);
B3 = (c^2 + px^2);
C3 = px*d3;
a = 1;
x_3 = circle(d2, A2, B2, C2, t_3 -4 -a, -a, x_2);
x_3(4) = d3;
for i = 1:length(t)
    current_t = t(i);
    if current_t < t_1
        x_star(:, i) = circle(d1, A1, B1, C1, current_t, 0, x0);
    elseif current_t < t_2
        x_star(:, i) = line(current_t, t_1, x_1);
    elseif current_t < t_3
        x_star(:, i) = circle(d2, A2, B2, C2, current_t - 4 -a, -a, x_2);
    else
        b = 2;
        x_star(:, i) = circle(d3, A3, B3, C3, current_t - t_3 -b, -b, x_3);
    end
end

end

