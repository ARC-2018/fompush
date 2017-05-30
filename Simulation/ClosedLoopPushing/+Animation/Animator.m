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
    
function [myMovie] = AnimateTrajectory(obj, frame_rate, file_name, objective_trajectory, simulated_trajectory, f_horizon_x_state, fh_considered_x, record)
    % Function to animate the pusher_slider simulation motion
    hFigure = figure('Color', 'w', 'OuterPosition', [0, 0, 960, 1080], ...
    'PaperPosition', [0, 0, 11, (6 / 8) * 11]);
    hold on 
    
%     frame1 = getframe(animation)
    if record == 1
        acc_factor = 1;
        speed_up = 1;
        %Create movie file
    else
%         acc_factor = 1;
        acc_factor = ceil((length(objective_trajectory) + 1) / 6);
        speed_up = 1;
    end
    videoname = strcat(file_name,'.avi');
    v = VideoWriter(videoname);
    v.FrameRate = frame_rate / acc_factor * speed_up;
    open(v);
    length(objective_trajectory)
    %Create label
    set(gcf, 'renderer', 'OpenGL');
    set(gca,'FontSize',20) % TODO: Idem here %Set size of axis font
    axis equal
    xlabel('x(m)', 'Interpreter', 'latex', 'FontSize', obj.font_size);
    ylabel('y(m)', 'Interpreter', 'latex', 'FontSize', obj.font_size);
    [x_lb_sim, x_ub_sim, y_lb_sim, y_ub_sim] = Models.QSPusherSlider.GetPlotLimits(simulated_trajectory);
    [x_lb_obj, x_ub_obj, y_lb_obj, y_ub_obj] = Models.QSPusherSlider.GetPlotLimits(objective_trajectory);
    xlim([min(x_lb_sim, x_lb_obj) max(x_ub_sim, x_ub_obj)]);
    ylim([min(y_lb_sim, y_lb_obj) max(y_ub_sim, y_ub_obj)]);
    %Go through mpc iterations
    for iteration = 1:acc_factor:length(objective_trajectory)
%     for iteration = [1, 600, 1200, 1800, 2400, length(objective_trajectory)]
        [sim_x_s, sim_y_s, sim_x_p, sim_y_p] = obj.data_provider.GetPusherSliderPolygons(simulated_trajectory(:, iteration));
        [obj_x_s, obj_y_s, obj_x_p, obj_y_p] = obj.data_provider.GetPusherSliderPolygons(objective_trajectory(:, iteration));
        if iteration == 1
            slider_line_width = 3.0;
            alpha = 1;
            MovingObjectiveSlider = patch(obj_x_s, obj_y_s, 'red', 'EdgeAlpha', 1, 'FaceAlpha', 1, 'EdgeColor', 'r', 'FaceColor', 'NONE', 'LineWidth', slider_line_width);
            MovingObjectivePusher = patch(obj_x_p, obj_y_p, 'red', 'EdgeAlpha', 1,'FaceAlpha', 1, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', [1,0,0] * 0.5, 'LineWidth', 0.1);
            MovingSimulatedSlider = patch(sim_x_s, sim_y_s, 'red', 'EdgeAlpha', 1, 'FaceAlpha', 1, 'EdgeColor', 'b', 'FaceColor', 'NONE', 'LineWidth', slider_line_width);
            MovingSimulatedPusher = patch(sim_x_p, sim_y_p, 'red', 'EdgeAlpha', alpha, 'FaceAlpha', alpha, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', [1,0,0] * 0.5, 'LineWidth', 0.1);
            MovingFiniteHorizon = cell(length(fh_considered_x{1}) + 1,1);
            for i = 1:length(fh_considered_x{iteration})
                if length(fh_considered_x{iteration}{i}) > 0
%                     MovingFiniteHorizon{i} = plot(fh_considered_x{iteration}{i}(1,:), fh_considered_x{iteration}{i}(2,:), 'Color', [0,0,0]);
                    MovingFiniteHorizon{i} = plot(fh_considered_x{iteration}{i}(1,1:3), fh_considered_x{iteration}{i}(2,1:3), 'Color', [0,0,0]);
                end
            end
            MovingFiniteHorizon{end} = plot(zeros(36, 1), zeros(36, 1), 'Color', [0,0,0,0.1]);
            plot(objective_trajectory(1, :), objective_trajectory(2, :), 'black');
            plot(simulated_trajectory(1, :), simulated_trajectory(2, :), 'blue');
            MovingChosenHorizon = plot(f_horizon_x_state{iteration}(1,:), f_horizon_x_state{iteration}(2,:), 'Color', [1,1,0], 'LineWidth', 2);
        else
            MovingObjectiveSlider.XData = obj_x_s;
            MovingObjectiveSlider.YData = obj_y_s;
            MovingObjectivePusher.XData = obj_x_p;
            MovingObjectivePusher.YData = obj_y_p;
            
            MovingSimulatedSlider.XData = sim_x_s;
            MovingSimulatedSlider.YData = sim_y_s;
            MovingSimulatedPusher.XData = sim_x_p;
            MovingSimulatedPusher.YData = sim_y_p;
            
            for i = 1:length(fh_considered_x{iteration})
                if length(fh_considered_x{iteration}{i}) > 0
%                     MovingFiniteHorizon{i}.XData = fh_considered_x{iteration}{i}(1,:);
%                     MovingFiniteHorizon{i}.YData = fh_considered_x{iteration}{i}(2,:);
                    MovingFiniteHorizon{i}.XData = fh_considered_x{iteration}{i}(1,1:3);
                    MovingFiniteHorizon{i}.YData = fh_considered_x{iteration}{i}(2,1:3);
                end
            end
%             MovingChosenHorizon.XData = f_horizon_x_state{iteration}(1,:);
%             MovingChosenHorizon.YData = f_horizon_x_state{iteration}(2,:);
            MovingChosenHorizon.XData = f_horizon_x_state{iteration}(1,1:3);
            MovingChosenHorizon.YData = f_horizon_x_state{iteration}(2,1:3);
            slider_line_width = 0.1;
            alpha = .4;
        end
        if ~record
            FiniteHorizon = cell(length(fh_considered_x{1}) + 1,1);
            for i = 1:length(fh_considered_x{iteration})
                if length(fh_considered_x{iteration}{i}) > 0
%                     FiniteHorizon{i} = plot(fh_considered_x{iteration}{i}(1,:), fh_considered_x{iteration}{i}(2,:), 'Color', [0,0,0]);
                    FiniteHorizon{i} = plot(fh_considered_x{iteration}{i}(1,1:3), fh_considered_x{iteration}{i}(2,1:3), 'Color', [0,0,0]);
                end
            end
%             ChozenHorizon = plot(f_horizon_x_state{iteration}(1,:), f_horizon_x_state{iteration}(2,:), 'Color', [1,1,0], 'LineWidth', 2);
            ChozenHorizon = plot(f_horizon_x_state{iteration}(1,1:3), f_horizon_x_state{iteration}(2,1:3), 'Color', [1,1,0], 'LineWidth', 2);
        end
        SimulatedSlider = patch(sim_x_s, sim_y_s, 'red', 'EdgeAlpha', .4, 'FaceAlpha', .4, 'EdgeColor', 'b', 'FaceColor', 'NONE', 'LineWidth', 0.1);
        SimulatedPusher = patch(sim_x_p, sim_y_p, 'red', 'EdgeAlpha', .4, 'FaceAlpha', .4, 'EdgeColor', [0,0,1] * 0.3, 'FaceColor', [1,0,0] * 0.5, 'LineWidth', 0.1);
        thisFrame = getframe(hFigure);
%         while true
%             w = waitforbuttonpress;
%             if w == 0
%                 break
%             end
%         end
        if record
            title(strcat('Velocity x', int2str(speed_up)))
            writeVideo(v,thisFrame);
        end
    end      
    close(v);     
%     close(hFigure)
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

