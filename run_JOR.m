% This is the code for the paper:
% Jointly optimized regressors for image super-resolution, 
% Dengxin Dai, Radu Timofte, Luc Van Gool, Eurographics 2015 
% 
% The code is built on top of the code developed by Radu Timofte and colleages for
% their ICCV13 work ANR and ACCV14 work A+. We also incorporated SRCNN into this code
% for comparison. All these methods can be compared by running this code directly  
% 
% search with <JOR> for places pertaining to JOR 
% make sure that you have the corresponding trained model: either download it at http://people.ee.ethz.ch/~daid/JOR/
% or train it with run_JOR_training.m
%
% Jan 10, 2015. Dengxin Dai, dai@vision.ee.ethz.ch, ddx2004@gmail.com

%%
addpath(genpath(pwd));  % the upscaling methods
addpath('vlfeat-0.9.20/mexmaci64');

imgscale = 1; % the scale reference we work with
flag = 1;       % flag = 0 - only bicubic methods, GR, ANR, JOR (Ours here), the other get the bicubic result by default
                % flag = 1 - all the methods are applied
upscaling = 3;  % {3, 4}, the magnification factor x3 or x4
input_dir = 'Set5';  %{Set5, Set14, B100, 'Tex136'}

if any(strcmp(input_dir, {'Set5', 'Set14'}))
    pattern = '*.bmp'; % Pattern to process
else 
    pattern = '*.jpg'; 
end
 
%% settings for alternatives such as ANR and A+ ...  

dict_sizes = [2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536];
neighbors = [1:1:12, 16:4:32, 40:8:64, 80:16:128, 256, 512, 1024];
setNum = [32 64 128 256 512]


clusterszA = 2048; % neighborhood size for A+

disp('The experiment corresponds to the results from Table 2 in the referenced [1] and [2] papers.');

disp(['The experiment uses ' input_dir ' dataset and aims at a magnification of factor x' num2str(upscaling) '.']);
if flag==1
    disp('All methods are employed : Bicubic, Yang et al., Zeyde et al., GR, ANR, NE+LS, NE+NNLS, NE+LLE, A+, SRCNN, JOR.');    
else
    disp('We run only for Bicubic, GR, ANR, A+, SRCNN and JOR, the other get the Bicubic result by default.');
end

fprintf('\n\n');

%%  super-resolution, note that JOR doesn't use dictionary 

window_num = [125,250,500,1000];
for w = 1:4 %window num to test 
    d=10    %1024
    conf_set = cell(1,10);    
    for s = 1:10 %use 10 dataset   
        s1 = 'img_';
        s2 = num2str(s);
        s3 = '_*.jpg';
        pattern_str = strcat(s1,s2,s3);
        conf.setNum = 10;
        conf.setSize = 100;
        conf.patch_size = [99,99];
        conf.scale = upscaling; % scale-up factor
        conf.level = 1; % # of scale-ups to perform
        conf.window = [3 3]; % low-res. window size
        conf.border = [1 1]; % border of the image (to ignore)
        conf.kmeans = 0;
        conf.kmeans_window = 0;
        conf.ksvd = 1;
        conf.patch_num = 11052;
        conf.window_num = window_num(1,w);
        conf.cluster_size = 5; % the cluster size used for k-means

        % High-pass filters for feature extraction (defined for upsampled low-res.)
        conf.upsample_factor = upscaling; % upsample low-res. into mid-res.
        O = zeros(1, conf.upsample_factor-1);
        G = [1 O -1]; % Gradient
        L = [1 O -2 O 1]/2; % Laplacian
        conf.filters = {G, G.', L, L.'}; % 2D versions
        conf.interpolate_kernel = 'bicubic';

        conf.overlap = [1 1]; % partial overlap (for faster training)
        if upscaling <= 2
            conf.overlap = [1 1]; % partial overlap (for faster training)
        end

        startt = tic;

%         conf.patch_num = size(winows,2);
        conf.overlap = conf.window - [1 1]; % full overlap scheme (for better reconstruction)    
        conf.trainingtime = toc(startt);
        toc(startt)            
        if dict_sizes(d) < 1024
            lambda = 0.01;
        elseif dict_sizes(d) < 2048
            lambda = 0.1;
        elseif dict_sizes(d) < 8192
            lambda = 1;
        else
            lambda = 5;
        end
        conf.filenames = glob(input_dir, pattern); % Cell array      

    %     conf.desc = {'Original', 'Bicubic', 'Yang et al.', ...
    %         'Zeyde et al.', 'GR', 'ANR', ...
    %         'NE+LS','NE+NNLS','NE+LLE', 'A+', 'SRCNN', 'JOR'};
        tag = [input_dir '_x' num2str(upscaling) '_' num2str(dict_sizes(d)) 'atoms'];
        conf.desc = {'Original','Bicubic','NE+LLE'};
        conf.results = {};
        windows = load_images_patches(glob('training_set', pattern_str),conf);%conf for different set
        conf = learn_dict(conf, windows(:,1:window_num(w)), dict_sizes(d));
        conf.points = [1:1:size(conf.dict_lores,2)];
        conf.pointslo = conf.dict_lores(:,conf.points);
        conf.pointsloPCA = conf.pointslo'*conf.V_pca';
        conf.result_dirImages = qmkdir([input_dir '/results_' tag]);
        conf.result_dirImagesRGB = qmkdir([input_dir '/results_' tag 'RGB']);
        conf.result_dir = qmkdir(['Results-' datestr(now, 'YYYY-mm-dd_HH-MM-SS')]);
        conf.result_dirRGB = qmkdir(['ResultsRGB-' datestr(now, 'YYYY-mm-dd_HH-MM-SS')]);
        
        disp(['Upscaling x' num2str(upscaling) ' ' input_dir ' with Zeyde dictionary of size = ' num2str(dict_sizes(d))]);

        mat_file = ['conf_Zeyde_' num2str(dict_sizes(d)) '_finalx' num2str(upscaling)];    

        t = cputime;    

        conf.countedtime = zeros(numel(conf.desc),numel(conf.filenames));

        res =[];
        for i = 1:numel(conf.filenames)
            f = conf.filenames{i};
            [p, n, x] = fileparts(f);
            [img, imgCB, imgCR] = load_images({f}); 
            if imgscale<1
                img = resize(img, imgscale, conf.interpolate_kernel);
                imgCB = resize(imgCB, imgscale, conf.interpolate_kernel);
                imgCR = resize(imgCR, imgscale, conf.interpolate_kernel);
            end
            sz = size(img{1});

            fprintf('%d/%d\t"%s" [%d x %d]\n', i, numel(conf.filenames), f, sz(1), sz(2));

            img = modcrop(img, conf.scale^conf.level);
            imgCB = modcrop(imgCB, conf.scale^conf.level);
            imgCR = modcrop(imgCR, conf.scale^conf.level);

                low = resize(img, 1/conf.scale^conf.level, conf.interpolate_kernel);
                if ~isempty(imgCB{1})
                    lowCB = resize(imgCB, 1/conf.scale^conf.level, conf.interpolate_kernel);
                    lowCR = resize(imgCR, 1/conf.scale^conf.level, conf.interpolate_kernel);
                end

            interpolated = resize(low, conf.scale^conf.level, conf.interpolate_kernel);
            if ~isempty(imgCB{1})
                interpolatedCB = resize(lowCB, conf.scale^conf.level, conf.interpolate_kernel);    
                interpolatedCR = resize(lowCR, conf.scale^conf.level, conf.interpolate_kernel);    
            end

            res{1} = interpolated;

            if flag == 1
                startt = tic;
                if 24 < dict_sizes(d)
                    res{2} = scaleup_NE_LLE(conf, low, 24);
                else
                    res{2} = scaleup_NE_LLE(conf, low, dict_sizes(d));
                end
                toc(startt)
                conf.countedtime(2,i) = toc(startt);    
            else
                res{2} = interpolated;
            end
            result = cat(3, img{1}, interpolated{1}, res{2}{1});

            result = shave(uint8(result * 255), conf.border * conf.scale);

            if ~isempty(imgCB{1})
                resultCB = interpolatedCB{1};
                resultCR = interpolatedCR{1};           
                resultCB = shave(uint8(resultCB * 255), conf.border * conf.scale);
                resultCR = shave(uint8(resultCR * 255), conf.border * conf.scale);
            end

            conf.results{i} = {};
            for j = 1:numel(conf.desc)            
                conf.results{i}{j} = fullfile(conf.result_dirImages, [n sprintf('[%d-%s]', j, conf.desc{j}) x]);            
                imwrite(result(:, :, j), conf.results{i}{j});

                conf.resultsRGB{i}{j} = fullfile(conf.result_dirImagesRGB, [n sprintf('[%d-%s]', j, conf.desc{j}) x]);
                if ~isempty(imgCB{1})
                    rgbImg = cat(3,result(:,:,j),resultCB,resultCR);
                    rgbImg = ycbcr2rgb(rgbImg);
                else
                    rgbImg = cat(3,result(:,:,j),result(:,:,j),result(:,:,j));
                end

                imwrite(rgbImg, conf.resultsRGB{i}{j});
            end        
            conf.filenames{i} = f;
        end   
        conf.duration = cputime - t;

    % Test performance
    scores = run_comparison(conf);
    conf.scores = scores;
%     process_scores_Tex(conf, scores,length(conf.filenames));
%      
%     run_comparisonRGB(conf); % provides color images and HTML summary
    conf_set{s} = conf;
    %    
    save([tag '_' mat_file '_results_imgscale_' num2str(imgscale)],'conf','scores');
    end
%     select_patches(conf);
    show_result(conf_set);
end
%
