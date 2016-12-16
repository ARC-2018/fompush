Foldername = 'Results/Experiments_2016_11_23/Videos/';
files = dir(strcat(Foldername, '*.MP4'));
Filenames = cell(length(files), 1);
for i = 1:length(files)
    Filenames{i} = files(i).name();
    n = length(Filenames{i});
    Filenames{i} = Filenames{i}(1:(n-4)); %crop .MP4 out
end
offsets_in_seconds = [-3.3 -2.3 -3.5 -3.1 -3 -3.5]; %TODO: Add a vector with all the offsets
run('Setup.m');
run('TestMatlabCalibration.m');
% for i = 1:1
for i = 1:length(Filenames)
    experiment_name = Filenames{i};
    offset_in_seconds = offsets_in_seconds(i);
 
    Data1 = load(strcat('Results/Experiments_2016_11_23/Data/', experiment_name, '.mat'));
    Data = Data1.Data;
    %% Initialize video
    video = VideoReader(strcat(Foldername, experiment_name,'.MP4'));
    frame_rate = video.FrameRate;
    %% Allocate space for the new video
    video_width = video.Width;
    video_height = video.Height;
    scaling_factor = video_width / length(imOrig);
    new_M = M * scaling_factor;
    new_M(3,3) = 1;
    % The trajectory in pixels is obtained, drawn and displayed
    new_video = struct('cdata',zeros(video_height, video_width ,3,'uint8'), 'colormap',[]);
    %% Read and modify each of the frames. Store the modified frame in a new video file
    k = 1;
    video_time = offset_in_seconds;
    N = length(Data.xMPC);
    while hasFrame(video)
        im = undistortImage(readFrame(video), cameraParameters);
        trajectory_index = find(Data.t > video_time, 1, 'first'); % The n*log(n) can be reduced to n, but it should go fast enough
        if video_time > Data.t(length(Data.t))
            trajectory_index = N;
        end
        x = 1000 * Data.xMPC{trajectory_index};
        y = 1000 * Data.yMPC{trajectory_index};
        trajectory = [x.' y.' ones(length(x), 1)];
        world_trajectory = RobotToWorldTrajectory(trajectory, 323.4, 100.2);
        camera_trajectory = WorldToCameraTrajectory(world_trajectory, new_M, planarR, translationVector);
        straight_line_trajectory = [200 0 1; 545 0 1];
        straight_line_world_trajectory = RobotToWorldTrajectory(straight_line_trajectory, 323.4, 100.2);
        straight_line_camera_trajectory = WorldToCameraTrajectory(straight_line_world_trajectory, new_M, planarR, translationVector);
        im = DrawLineBetweenTwoPoints(straight_line_camera_trajectory, im);
        new_video(k).cdata = DrawTrajectoryOnImage(camera_trajectory, im);
        k = k+1;
        video_time = video_time + 1/frame_rate;
    end
    video_writer = VideoWriter(strcat('Results/Experiments_2016_11_23/Videos/', experiment_name));
    video_writer.FrameRate = video.FrameRate;
    open(video_writer);
    writeVideo(video_writer, new_video);
    close(video_writer);
end