% This is the training code of JOR for the paper:
% Jointly optimized regressors for image super-resolution, 
% Dengxin Dai, Radu Timofte, Luc Van Gool, Eurographics 2015 

% written by Dengxin Dai at ETH Zurich
% dai@vision.ee.ethz.ch, ddx2004@gmail.com

% dependency: vl_feat library

% parameters you might want to tweak: 
%     upscaling: upscaling factor 
%     tr_num   : number of LR and HR samples 
%     scale_num: number of scales for downsampling to collect samples 
%     downsample_factor: downsampling factor
%     em_inum: number of iterations for EM

%% 
clear all;  
addpath(genpath(pwd));
                % flag = 1 - all the methods are applied
upscaling =3;   % the magnification factor x2, x3, x4...
win_size = 3; 

% path to the training images
input_dir = './libs/AplusCodes_SR/CVPR08-SR/Data/Training';
pattern = '*.bmp'; % Pattern to process

dict_sizes = [16 32 64 128 ];  % 

for d=2  
      
    disp(['Training dictionary of size ' num2str(dict_sizes(d)) ' using EM approach...']);
        % Simulation settings
    JOR.scale = upscaling; % scale-up factor
    JOR.level = 1; % # of scale-ups to perform
    JOR.window = [win_size win_size]; % low-res. window size
    JOR.border = [1 1]; % border of the image (to ignore)

    % High-pass filters for feature extraction (defined for upsampled low-res.)
    JOR.upsample_factor = upscaling; % upsample low-res. into mid-res.
    O = zeros(1, JOR.upsample_factor-1);
    G = [ 1 O -1 ]; % Gradient
    L = [1 O -2 O 1]/2; % Laplacian
    JOR.filters = {G, G.', L, L.'}; % 2D versions
    JOR.interpolate_kernel = 'bicubic';

    JOR.overlap = [2 2]; % partial overlap (for faster training)
    if upscaling <= 2
        JOR.overlap = [2 2]; % partial overlap (for faster training)
    end

    JOR.overlap = JOR.window - [1 1]; % full overlap scheme (for better reconstruction)    
         
    JOR.filenames = glob(input_dir, pattern); % Cell array  
 

  %% collecting samples   
    vl_setup; 
    dict_size = dict_sizes(d);
  
    tr_num = 5000000;  % number of samples 
    % tweak these two parameters to get enough samples 
    scale_num = 20; 
    downsample_factor = 0.98;
    [plores phires JOR.V_pca] = collectSamplesMulScales(JOR, load_images(glob(input_dir, ['*', pattern])), scale_num, downsample_factor);

    if size(plores,2) > tr_num     
        s_factor = size(plores,2)/tr_num;
        plores = plores(:,1:s_factor:end);
        phires = phires(:,1:s_factor:end);
    end
    number_samples = size(plores,2);
    
    % l2 normalize LR patches, and scale the corresponding HR patches
    l2 = sum(plores.^2).^0.5+eps;
    l2n = repmat(l2,size(plores,1),1);   
%     l2(l2<0.1) = 1;
    plores = plores./l2n;
    phires = phires./repmat(l2,size(phires,1),1);
    clear l2
    clear l2n
    JOR.lowpatches = plores; 
    JOR.highpatches = phires; 
    clear plores
    clear phires
    
    %% JOR training 
    em_inum = 20;    % # of iterations of EM
    [JOR.PPs JOR.labels, JOR.error_table] = Learn_Projections_EM(JOR.lowpatches, JOR.highpatches, dict_size, em_inum);

%% clean and save the model
    JOR.upscaling = JOR.scale;
    JOR = rmfield(JOR, 'scale');
    JOR = rmfield(JOR, 'filenames');
    JOR = rmfield(JOR, 'highpatches');
    
    mat_file = ['JOR_model_' num2str(dict_sizes(d)) '_x' num2str(upscaling) '_' num2str(tr_num) '.mat'];    
    save('-v7.3', ['./models/' mat_file], 'JOR');                       
    
    disp('training finished !');
end
