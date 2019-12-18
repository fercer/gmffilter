% Test minimizing the MSE optimizing T, L, and sigma
clear all
close all
clc

%% Set the environment variables:
%
DB_dir = 'E:/test_data/Angios_134/';
Mask_tools_dir = 'E:/Apps/detection_tools/matlab_version/mask_tools';
Performance_tools_dir = 'E:/Apps/detection_tools/matlab_version/performance_tools';
%}
%{
DB_dir = 'C:/test_data/';
Mask_tools_dir = 'C:/Users/cimat/detection_tools/matlab_version/mask_tools';
Performance_tools_dir = 'C:/Users/cimat/detection_tools/matlab_version/performance_tools';
%}
addpath(DB_dir);
addpath(Mask_tools_dir);
addpath(Performance_tools_dir);
addpath('.\Image_tools\');


%% Load images
global training_set_size testing_set_size training_set testing_set 
global mask img_gt raw_img fft_img
global cols rows old_rows old_cols offset_x offset_y
global K beta delta

training_set_size = 100;
testing_set_size = 30;

total_images_size = training_set_size + testing_set_size;

mask = zeros(300, 300, total_images_size);
img_gt = zeros(300, 300, total_images_size);
y_target = zeros(300, 300, total_images_size);
raw_img = zeros(300, 300, total_images_size);
fft_img = zeros(512, 512, total_images_size);
for img_i = 1:total_images_size
    img_filename = sprintf('%d.png',img_i);
    img = double(imread(img_filename)) / 255.0;
    mask_idx = compute_FOV_mask(img);
    mask(:,:,img_i) = mask_corners(mask_idx);
    img = fill_black_area(img, mask_idx);
    raw_img(:,:,img_i) = img;

    img_gt_filename = sprintf('%d_gt.png',img_i);
    img_gt(:,:,img_i) = double(imread(img_gt_filename)) / 255.0;
    [~, y_target_i] = image_preprocessing(img);
    y_target(:,:,img_i) = y_target_i;
    [old_rows, old_cols] = size(img);

    % Pad the image to have sizes of 2-powered multiple
    img = padImage(raw_img(:,:,img_i));
    fft_img(:,:,img_i) = fft2(img);
end

[rows, cols] = size(img);
offset_y = rows/2 - old_rows/2;
if offset_y < 1
    offset_y = 1;
end

offset_x = cols/2 - old_cols/2;
if offset_x  < 1
    offset_x = 1;
end

%% Separate the training, validation and testing sets
% From the training set
validation_set_percentage = 0.0;

% Separate the images in each set
training_set = 1:training_set_size;
testing_set = (1:testing_set_size) + training_set_size;

% Extract the validation set from the training set:
validation_set_size = floor(validation_set_percentage * training_set_size);
validation_set = randsample(training_set, validation_set_size, false);
training_set = setdiff(training_set, validation_set);
training_set_size = training_set_size - validation_set_size;


%% Set the GMF initial parameters:
% For the step functions:
delta = 1e-0;

% For the softmax functions:
beta = 1e3;

%% Set the GMF global variables
K = 12;
x = [2.82; 13; 15; 1; 0];

GMF_resp_ini = gmf_fun_resp(x);
dGMF_resp_ini = dgmf_fun_resp(x);
GMF_resp_ini_tst = gmf_fun_resp_tst(x);

%% Test the LM algorithm
x_ini = [5.0; 15.0; 15.0; 1; 0];
x_min = [0.1; 0.5; 1; 0.001; -1.0];
x_max = [10; 30; 512; 1000; 1.0];

x_opt = levenberg_marquardt(x_ini, reshape(img_gt(:,:,training_set), [training_set_size*old_rows*old_cols, 1]), 'gmf_fun_resp', 'dgmf_fun_resp', x_min, x_max);
x_opt

%x_opt = [3.1803;22.8712;15.0695;25.2787;-0.0194];
%% Test the LM result:
GMF_resp_opt = gmf_fun_resp(x_opt);
GMF_resp_opt_tst = gmf_fun_resp_tst(x_opt);

figure
subplot(3, 2, 1)
imshow(reshape(GMF_resp_ini((1:(old_cols*old_rows))), [old_rows, old_cols]),[])

subplot(3, 2, 3)
imshow(reshape(GMF_resp_ini((1:(old_cols*old_rows))+old_cols*old_rows), [old_rows, old_cols]),[])

subplot(3, 2, 5)
imshow(reshape(GMF_resp_ini((1:(old_cols*old_rows))+2*old_cols*old_rows), [old_rows, old_cols]),[])

subplot(3, 2, 2)
imshow(reshape(GMF_resp_opt((1:(old_cols*old_rows))), [old_rows, old_cols]),[])

subplot(3, 2, 4)
imshow(reshape(GMF_resp_opt((1:(old_cols*old_rows))+old_cols*old_rows), [old_rows, old_cols]),[])

subplot(3, 2, 6)
imshow(reshape(GMF_resp_opt((1:(old_cols*old_rows))+2*old_cols*old_rows), [old_rows, old_cols]),[])

[~, ~, Az_ini_trn] = run_ROC(mask(:,:,training_set), GMF_resp_ini, img_gt(:,:,training_set))
[~, ~, Az_opt_trn] = run_ROC(mask(:,:,training_set), GMF_resp_opt, img_gt(:,:,training_set))

[~, ~, Az_ini_tst] = run_ROC(mask(:,:,testing_set), GMF_resp_ini_tst, img_gt(:,:,testing_set))
[~, ~, Az_opt_tst] = run_ROC(mask(:,:,testing_set), GMF_resp_opt_tst, img_gt(:,:,testing_set))


%% Test the LM for cross entropy minimization algorithm:
x_opt = levenberg_marquardt_cross_entropy(x, reshape(img_gt(:,:,training_set), [old_cols*old_rows*training_set_size, 1]), 'gmf_fun_resp', 'dgmf_fun_resp', x_min, x_max);
x_opt

%% Test the LM/CE result:
fx_opt = gmf_fun_resp(x_opt);
GMF_resp_opt = GMF_resp;

figure
subplot(3, 2, 1)
imshow(GMF_resp_ini(:,:,2),[])

subplot(3, 2, 3)
imshow(GMF_resp_ini(:,:,4),[])

subplot(3, 2, 5)
imshow(GMF_resp_ini(:,:,8),[])

subplot(3, 2, 2)
imshow(GMF_resp_ini(:,:,2),[])

subplot(3, 2, 4)
imshow(GMF_resp_opt(:,:,4),[])

subplot(3, 2, 6)
imshow(GMF_resp_opt(:,:,8),[])

%[~, ~, Az_ini] = run_ROC(mask(:,:,training_set), GMF_resp_ini,img_gt(:,:,training_set))
%[~, ~, Az_opt] = run_ROC(mask(:,:,training_set), GMF_resp_opt,img_gt(:,:,training_set))
