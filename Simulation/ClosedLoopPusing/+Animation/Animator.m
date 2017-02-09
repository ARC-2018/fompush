classdef Animator
%ANIMATOR Class to set parameters for the animation. It also includes all
%the possible animation functions.

properties (Access = private)
    font_size  = 25;
    line_size  = 15;
    line_width = 2;
    data_provider;
end 

methods
function obj = Animator(data_provider)
    obj.data_provider = data_provider;
end
    
function [] = AnimateTrajectory(obj, frame_rate, file_name, objective_trajectory, simulated_trajectory)
    % Function to animate the pusher_slider simulation motion
    animation = figure('Color', 'w', 'OuterPosition', [0, 0, 960, 1080], ...
    'PaperPosition', [0, 0, 11, (6 / 8) * 11]);
    acc_factor = 5;
    set(gcf,'Renderer','OpenGL'); % TODO: Check if it actually does anything
    set(gca,'FontSize',20) % TODO: Idem here %Set size of axis font
    axis equal
    %Create label
    xlabel('x(m)', 'Interpreter', 'latex', 'FontSize', obj.font_size);
    ylabel('y(m)', 'Interpreter', 'latex', 'FontSize', obj.font_size);
    % zlabel('z(m)', 'Interpreter', 'latex', 'FontSize', 16);
    % TODO: Get xlim and ylim from getter functions
    %Create movie file
    videoname = strcat(file_name,'.avi');
    v = VideoWriter(videoname);
    v.FrameRate = frame_rate;
    open(v);
    [x_lb_sim, x_ub_sim, y_lb_sim, y_ub_sim] = Models.QSPusherSlider.GetPlotLimits(simulated_trajectory);
    [x_lb_obj, x_ub_obj, y_lb_obj, y_ub_obj] = Models.QSPusherSlider.GetPlotLimits(objective_trajectory);
    xlim([min(x_lb_sim, x_lb_obj) max(x_ub_sim, x_ub_obj)]);
    ylim([min(y_lb_sim, y_lb_obj) max(y_ub_sim, y_ub_obj)]);
    %Go through mpc iterations
    for iteration = 1:acc_factor:length(objective_trajectory)
        [sim_x_s, sim_y_s, sim_x_p, sim_y_p] = obj.data_provider.GetPusherSliderPolygons(simulated_trajectory(:, iteration));
        [obj_x_s, obj_y_s, obj_x_p, obj_y_p] = obj.data_provider.GetPusherSliderPolygons(objective_trajectory(:, iteration));
        if iteration == 1
            ObjectiveSlider = patch(obj_x_s, obj_y_s, 'red', 'EdgeAlpha', 1, 'FaceAlpha', 1, 'EdgeColor', 'r', 'FaceColor', 'NONE', 'LineWidth', 0.1);
            hold on 
            ObjectivePusher = patch(obj_x_p, obj_y_p, 'red', 'EdgeAlpha', 1,'FaceAlpha', 1, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', [1,0,0] * 0.5, 'LineWidth', 0.1);
            slider_line_width = 3.0;
            alpha = 1;
        else
            ObjectiveSlider.XData = obj_x_s;
            ObjectiveSlider.YData = obj_y_s;
            ObjectivePusher.XData = obj_x_p;
            ObjectivePusher.YData = obj_y_p;
            slider_line_width = 0.1;
            alpha = .2;
        end
        SimulatedSlider = patch(sim_x_s, sim_y_s, 'red', 'EdgeAlpha', alpha, 'FaceAlpha', alpha, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', 'NONE', 'LineWidth', slider_line_width);
        hold on 
        SimulatedPusher = patch(sim_x_p, sim_y_p, 'red', 'EdgeAlpha', alpha, 'FaceAlpha', alpha, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', [1,0,0] * 0.5, 'LineWidth', 0.1);
        frame = getframe(animation);
        writeVideo(v,frame);
    end           
    close(v);
end

function [] = AnimateTracking(obj, file_name, simulated_trajectory, objective_points)
    % Function to animate the pusher_slider simulation motion
    animation = figure('Color', 'w', 'OuterPosition', [0, 0, 960, 1080], ...
    'PaperPosition', [0, 0, 11, (6 / 8) * 11]);
    acc_factor = 5;
    set(gcf,'Renderer','OpenGL'); % TODO: Check if it actually does anything
    set(gca,'FontSize',20) % TODO: Idem here %Set size of axis font
    axis equal
    %Create label
    xlabel('x(m)', 'Interpreter', 'latex', 'FontSize', obj.font_size);
    ylabel('y(m)', 'Interpreter', 'latex', 'FontSize', obj.font_size);
    % TODO: Get xlim and ylim from getter functions
    [x_lb, x_up, y_lb, y_ub] = Models.QSPusherSlider.GetPlotLimits(simulated_trajectory);
    xlim([x_lb x_up]);
    ylim([y_lb y_ub]);
    for iteration = 1:acc_factor:length(simulated_trajectory)
        [sim_x_s, sim_y_s, sim_x_p, sim_y_p] = obj.data_provider.GetPusherSliderPolygons(simulated_trajectory(:, iteration));
        if iteration == 1
            slider_line_width = 3.0;
            alpha = 1;
        else
            slider_line_width = 0.1;
            alpha = .2;
        end
        Slider = patch(sim_x_s, sim_y_s, 'red', 'EdgeAlpha', alpha, 'FaceAlpha', alpha, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', 'NONE', 'LineWidth', slider_line_width);
        hold on 
        Pusher = patch(sim_x_p, sim_y_p, 'red', 'EdgeAlpha', alpha, 'FaceAlpha', alpha, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', [1,0,0] * 0.5, 'LineWidth', 0.1);
    end
    h = scatter(objective_points(1,:), objective_points(2, :), 100, 'b', '^', 'filled');
    legend(h,'Target');
    saveas(animation, file_name, 'epsc');
end

end
    
end

