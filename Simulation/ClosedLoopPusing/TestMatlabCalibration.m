%% Initialize the trajectory. Hopefully, only these lines will be changed
N = 200; % Select number of points in the trajectory
t = linspace(0, 2*pi, N);
% x = 16 * sin(t).^3; % Store a vector for x and one for y with the positions
% y = 13 * cos(t) - 5 * cos(2*t) - 2 * cos(3*t) - cos(4*t); % The values for
x = 40 * t;
y = 0;
% X and Y MUST be provided in mm. The size of one square is currently 26.5

%% Image calibration setup
numImages = 13; % Number of images to calibrate
files = cell(1, numImages);
for i = 1:numImages % The images are stored in a cell
   files{i} = fullfile(sprintf('img/image%d.jpg', i));
end
squareSize = 26.5; % in millimeters
%% Image and world points generation
% First the points are detected in the image. The squares per side are also
% computed
[imagePoints, boardSize] = detectCheckerboardPoints(files);
% From the number of squares and their size, the world points are computed
% The reference frame starts in one of the black outer vertices
worldPoints = generateCheckerboardPoints(boardSize, squareSize);
%% Intrinsics computation
cameraParams = estimateCameraParameters(imagePoints, worldPoints);
%% Extrinsics computation for the recording position
% The image taken from the recording position is loaded and displayed
imOrig = imread(fullfile('img/image14.jpg'));
im = undistortImage(imOrig, cameraParams); % The image is undistorted
% The procedure was repeated in the Matlab website, presumably to allow the
% usage of different checkerboards for intrinsic and extrinsic calibration
% or if the intrinsics are already known. We can ignore the previous lines
% if we use the ros calibration script for intrinsics, or we can remove
% this line if we use Matlab and the same board for intrinsic and extrinsic
[imagePoints, boardSize] = detectCheckerboardPoints(im);
% R and T are computed
[rotationMatrix, translationVector] = extrinsics(imagePoints, worldPoints, cameraParams);
translationVector = translationVector.'; % My code used vertical vectors
M = cameraParams.IntrinsicMatrix.'; % Matlab returns a transposed M
%% Trajectory transformation
trajectory = zeros(size(worldPoints,1), 3);
trajectory(:,1:2) = worldPoints; % The trajectory points are initialized in the
% trajectory(:,2) = y.'; % format used in my functions
trajectory(:,3) = ones(size(worldPoints,1),1);
new_theta = 0.175;
rotationMatrix = [cos(new_theta), -sin(new_theta) 0; sin(new_theta), cos(new_theta) 0; 0 0 1] * rotationMatrix;
planarR = rotationMatrix(:,1:2); % Everything is on the same plane so we can set Z = 0
% The trajectory in pixels is obtained, drawn and displayed
camera_trajectory = WorldToCameraTrajectory(trajectory, M, planarR, translationVector);
new_I = DrawTrajectoryOnImage(camera_trajectory, im, [255 0 0]);
figure; imshow(new_I);