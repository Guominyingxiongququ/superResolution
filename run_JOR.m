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

figure_num = 1;
imgscale = 1; % the scale reference we work with
flag = 1;       % flag = 0 - only bicubic methods, GR, ANR, JOR (Ours here), the other get the bicubic result by default
                % flag = 1 - all the methods are applied
upscaling = 3;  % {3, 4}, the magnification factor x3 or x4
% input_dir = 'test_set';  %{Set5, Set14, B100, 'Tex136'}
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
test_filenames = glob(input_dir, pattern); % Cell array
window_num = [125,250,500,1000];
figure_num = 10;
for w = 1:4 %window num to test 
    d=10    %1024
    
    % insert other methods
    conf_set = cell(1,10); 
    %active learning iteration
    for iter = 1
        filenames = cell(1,10);
        for file = 1:10
%             filenames{file} = test_filenames{file+(iter-1)*10:iter*10};
            filenames{file} = test_filenames;
        end
        conf.filenames = filenames;
        
        % use which patch number of dataset from 1 to 10
        for s = 1 %use 10 dataset   
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
            conf.word_num = 1000; % set word number to 1000
            conf.cluster_size = 5; % the cluster size used for k-means

            % High-pass filters for feature extraction (defined for upsampled low-res.)
            conf.upsample_factor = upscaling; % upsample low-res. into mid-res.
            O = zeros(1, conf.upsample_factor-1);
            G = [1 O -1]; % Gradient
            L = [1 O -2 O 1]/2; % Laplacian
            conf.filters = {G, G.', L, L.'}; % 2D versions
            conf.interpolate_kernel = 'bicubic';
            conf.overlap = [1 1]; % partial overlap (for faster training)
            startt = tic;
    %  conf.patch_num = size(winows,2);
            
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
            
            windows = load_images_windows(glob('training_set', pattern_str),conf);%conf for different set
            conf = learn_dict(conf, windows(:,1:window_num(w)), dict_sizes(d));
            
            %% GR
            if dict_sizes(d) < 10000
                conf.ProjM = inv(conf.dict_lores'*conf.dict_lores+lambda*eye(size(conf.dict_lores,2)))*conf.dict_lores';    
                conf.PP = (1+lambda)*conf.dict_hires*conf.ProjM;
            else
                % here should be an approximation
                conf.PP = zeros(size(conf.dict_hires,1), size(conf.V_pca,2));
                conf.ProjM = [];
            end
            
            %%config and window
%             conf.filenames = glob(input_dir, pattern); % Cell array
            conf.desc = {'Original', 'Bicubic', 'Yang et al.', ...
                'Zeyde et al.', 'GR', 'ANR', ...
                'NE+LS','NE+NNLS','NE+LLE', 'A+', 'SRCNN' };
            conf.results = {};  
            tag = [input_dir '_x' num2str(upscaling) '_' num2str(dict_sizes(d)) 'atoms'];
            conf.points = [1:1:size(conf.dict_lores,2)];
            conf.pointslo = conf.dict_lores(:,conf.points);
            conf.pointsloPCA = conf.pointslo'*conf.V_pca';
            conf.result_dirImages = qmkdir([input_dir '/results_' tag]);
            conf.result_dirImagesRGB = qmkdir([input_dir '/results_' tag 'RGB']);
            conf.result_dir = qmkdir(['Results-' datestr(now, 'YYYY-mm-dd_HH-MM-SS')]);
            conf.result_dirRGB = qmkdir(['ResultsRGB-' datestr(now, 'YYYY-mm-dd_HH-MM-SS')]);
            conf.entropy = cell(1,5);
            disp(['Upscaling x' num2str(upscaling) ' ' input_dir ' with Zeyde dictionary of size = ' num2str(dict_sizes(d))]);

            mat_file = ['conf_Zeyde_' num2str(dict_sizes(d)) '_finalx' num2str(upscaling)];    
            
            %% precomputer for ANR
            % precompute for ANR the anchored neighborhoods and the projection matrices for
            % the dictionary 
            
            conf.PPs = [];    
            if  size(conf.dict_lores,2) < 40
                clustersz = size(conf.dict_lores,2);
            else
                clustersz = 40;
            end
            D = abs(conf.pointslo'*conf.dict_lores);    

            for i = 1:length(conf.points)
                [vals idx] = sort(D(i,:), 'descend');
                if (clustersz >= size(conf.dict_lores,2)/2)
                    conf.PPs{i} = conf.PP;
                else
                    Lo = conf.dict_lores(:, idx(1:clustersz));        
                    conf.PPs{i} = 1.01*conf.dict_hires(:,idx(1:clustersz))*inv(Lo'*Lo+0.01*eye(size(Lo,2)))*Lo';    
                end
            end    

            ANR_PPs = conf.PPs; % store the ANR regressors

            save([tag '_' mat_file '_ANR_projections_imgscale_' num2str(imgscale)],'conf');
            t = cputime;    
            
            %% A+ computing the regressors
            Aplus_PPs = [];
        
            fname = ['Aplus_x' num2str(upscaling) '_' num2str(dict_sizes(d)) 'atoms' num2str(clusterszA) 'nn_5mil.mat'];
    
%             if exist(fname,'file')
%                load(fname);
%             else
            disp('Compute A+ regressors');
            ttime = tic;
            tic
           [plores phires] = collectSamplesScales(conf, windows(:,1:window_num(w)), 12, 0.98);  

            if size(plores,2) > 5000000                
                plores = plores(:,1:5000000);
                phires = phires(:,1:5000000);
            end
            number_samples = size(plores,2);

            % l2 normalize LR patches, and scale the corresponding HR patches
            l2 = sum(plores.^2).^0.5+eps;
            l2n = repmat(l2,size(plores,1),1);    
            l2(l2<0.1) = 1;
            plores = plores./l2n;
            phires = phires./repmat(l2,size(phires,1),1);
            clear l2
            clear l2n

            llambda = 0.1;

            for i = 1:size(conf.dict_lores,2)
                D = pdist2(single(plores'),single(conf.dict_lores(:,i)'));
                [~, idx] = sort(D);                
                Lo = plores(:, idx(1:clusterszA));                                    
                Hi = phires(:, idx(1:clusterszA));
                Aplus_PPs{i} = Hi*inv(Lo'*Lo+llambda*eye(size(Lo,2)))*Lo'; 
            %Aplus_PPs{i} = Hi*(inv(Lo*Lo'+llambda*eye(size(Lo,1)))*Lo)'; 
            end        
            clear plores
            clear phires

            ttime = toc(ttime);        
            save(fname,'Aplus_PPs','ttime', 'number_samples');   
            toc
%             end    
            
            %% A+ (0.5mil) computing the regressors with 0.5 milion training samples
            Aplus05_PPs = [];    

            fname = ['Aplus_x' num2str(upscaling) '_' num2str(dict_sizes(d)) 'atoms' num2str(clusterszA) 'nn_05mil.mat'];    
    
%             if exist(fname,'file')
%                load(fname);
%             else
           disp('Compute A+ (0.5 mil) regressors');
           ttime = tic;
           tic
           [plores phires] = collectSamplesScales(conf, windows(:,1:window_num(w)), 1,1);  

            if size(plores,2) > 500000                
                plores = plores(:,1:500000);
                phires = phires(:,1:500000);
            end
            number_samples = size(plores,2);

            % l2 normalize LR patches, and scale the corresponding HR patches
            l2 = sum(plores.^2).^0.5+eps;
            l2n = repmat(l2,size(plores,1),1);      
            l2(l2<0.1) = 1;
            plores = plores./l2n;
            phires = phires./repmat(l2,size(phires,1),1);
            clear l2
            clear l2n

            llambda = 0.1;

            for i = 1:size(conf.dict_lores,2)
                D = pdist2(single(plores'),single(conf.dict_lores(:,i)'));
                [~, idx] = sort(D);                
                Lo = plores(:, idx(1:clusterszA));                                    
                Hi = phires(:, idx(1:clusterszA));
                Aplus05_PPs{i} = Hi*inv(Lo'*Lo+llambda*eye(size(Lo,2)))*Lo'; 
            end        
            clear plores
            clear phires

            ttime = toc(ttime);        
            save(fname,'Aplus05_PPs','ttime', 'number_samples');   
            toc
%             end            
    
            %% load the A+ (16 atoms) for comparison results
            conf16 = [];       
            fname = ['Aplus_x' num2str(upscaling) '_16atoms' num2str(clusterszA) 'nn_05mil.mat'];
            fnamec = ['Set14_x' num2str(upscaling) '_16atoms_conf_Zeyde_16_finalx' num2str(upscaling) '_ANR_projections_imgscale_' num2str(imgscale) '.mat']; 
%             if exist(fname,'file') && exist(fnamec,'file')
%                kk = load(fnamec);
%                conf16 = kk.conf;       
%                kk = load(fname);       
%                conf16.PPs = kk.Aplus05_PPs;
%                clear kk
%             end
            %%   build result folder 
            conf.result_dirImages = qmkdir([input_dir '/results_' tag]);
            conf.result_dirImagesRGB = qmkdir([input_dir '/results_' tag 'RGB']);
            conf.result_dir = qmkdir(['Results-' datestr(now, 'YYYY-mm-dd_HH-MM-SS')]);
            conf.result_dirRGB = qmkdir(['ResultsRGB-' datestr(now, 'YYYY-mm-dd_HH-MM-SS')]);

            t = cputime;    

            conf.countedtime = zeros(numel(conf.desc),numel(conf.filenames));
            %% prediction part
            res =[];
            filenames_list = conf.filenames{1};
            for i = 1:numel(conf.filenames{1})
                f = filenames_list{i};
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

                if (flag == 1) && (dict_sizes(d) == 1024) && (upscaling==3)
                    startt = tic;
                    res{2} = {yima(low{1}, upscaling)};                        
                    toc(startt)
                    conf.countedtime(2,i) = toc(startt);
                else
                    res{2} = interpolated;
                end

                if (flag == 1)
                    startt = tic;
                    res{3} = scaleup_Zeyde(conf, low);
                    toc(startt)
                    conf.countedtime(3,i) = toc(startt);    
                else
                    res{3} = interpolated;
                end

                %if flag == 1
                startt = tic;
                res{4} = scaleup_GR(conf, low);
                toc(startt)
                conf.countedtime(4,i) = toc(startt);    
                %else
                    %res{4} = interpolated;
                %end

                startt = tic;
                conf.PPs = ANR_PPs;
                res{5} = scaleup_ANR(conf, low);
                toc(startt)
                conf.countedtime(5,i) = toc(startt);    

                if flag == 1
                    startt = tic;
                    if 12 < dict_sizes(d)
                        res{6} = scaleup_NE_LS(conf, low, 12);
                    else
                        res{6} = scaleup_NE_LS(conf, low, dict_sizes(d));
                    end
                    toc(startt)
                    conf.countedtime(6,i) = toc(startt);    
                else
                    res{6} = interpolated;
                end

                if flag == 1
                    startt = tic;
                    if 24 < dict_sizes(d)
                        res{7} = scaleup_NE_NNLS(conf, low, 24);
                    else
                        res{7} = scaleup_NE_NNLS(conf, low, dict_sizes(d));
                    end
                    toc(startt)
                    conf.countedtime(7,i) = toc(startt);    
                else
                    res{7} = interpolated;
                end

                if flag == 1
                    startt = tic;
                    if 24 < dict_sizes(d)
                        res{8} = scaleup_NE_LLE(conf, low, 24);
                    else
                        res{8} = scaleup_NE_LLE(conf, low, dict_sizes(d));
                    end
                    toc(startt)
                    conf.countedtime(8,i) = toc(startt);    
                else
                    res{8} = interpolated;
                end
            
                % A+
                if ~isempty(Aplus_PPs)
                    fprintf('A+\n');
                    conf.PPs = Aplus_PPs;
                    startt = tic;
                    res{9} = scaleup_ANR(conf, low);
                    toc(startt)
                    conf.countedtime(9,i) = toc(startt);    
                else
                    res{9} = interpolated;
                end   

                % SRCNN
                startt = tic;
                res{10} = scaleup_SRCNN(conf, low);
                toc(startt)
                conf.countedtime(10,i) = toc(startt);   
            
            
            
%             conf.countedtime = zeros(numel(conf.desc),numel(conf.filenames));
% 
%             res =[];
%             for i = 1:numel(conf.filenames)
%                 f = conf.filenames{i};
%                 [p, n, x] = fileparts(f);
%                 [img, imgCB, imgCR] = load_images({f}); 
%                 if imgscale<1
%                     img = resize(img, imgscale, conf.interpolate_kernel);
%                     imgCB = resize(imgCB, imgscale, conf.interpolate_kernel);
%                     imgCR = resize(imgCR, imgscale, conf.interpolate_kernel);
%                 end
%                 sz = size(img{1});
% 
%                 fprintf('%d/%d\t"%s" [%d x %d]\n', i, numel(conf.filenames), f, sz(1), sz(2));
% 
%                 img = modcrop(img, conf.scale^conf.level);
%                 imgCB = modcrop(imgCB, conf.scale^conf.level);
%                 imgCR = modcrop(imgCR, conf.scale^conf.level);
% 
%                     low = resize(img, 1/conf.scale^conf.level, conf.interpolate_kernel);
%                     if ~isempty(imgCB{1})
%                         lowCB = resize(imgCB, 1/conf.scale^conf.level, conf.interpolate_kernel);
%                         lowCR = resize(imgCR, 1/conf.scale^conf.level, conf.interpolate_kernel);
%                     end
% 
%                 interpolated = resize(low, conf.scale^conf.level, conf.interpolate_kernel);
%                 if ~isempty(imgCB{1})
%                     interpolatedCB = resize(lowCB, conf.scale^conf.level, conf.interpolate_kernel);    
%                     interpolatedCR = resize(lowCR, conf.scale^conf.level, conf.interpolate_kernel);    
%                 end
% 
%                 res{1} = interpolated;
% 
%                 if flag == 1
%                     startt = tic;
%                     if 24 < dict_sizes(d)
%                         res{2} = scaleup_NE_LLE(conf, low, 24);
%                     else
%                         res{2} = scaleup_NE_LLE(conf, low, dict_sizes(d));
%                     end
%                     toc(startt)
%                     conf.countedtime(2,i) = toc(startt);    
%                 else
%                     res{2} = interpolated;
%                 end
                
                 result = cat(3, img{1}, interpolated{1}, res{2}{1}, res{3}{1}, ...
                    res{4}{1}, res{5}{1}, res{6}{1}, res{7}{1}, res{8}{1}, ...
                    res{9}{1}, res{10}{1});

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
        
        %% show entropy plot
        patch_list = select_patches(conf);
        conf.window = [3,3];
        show_patches(patch_list,conf,100,3,figure_num);
        figure_num = figure_num + 1;
        conf.window = [99,99];
        entropy = zeros(1,5);
        PSNR_score = zeros(1,5);
        for ranking = 1:5 
            [window_list,score_list] = select_windows(conf,ranking);
            entropy(ranking) = get_entropy(window_list);
            PSNR_score(ranking) = mean(score_list);
            show_patches(window_list,conf,10,1,figure_num);
            figure_num = figure_num + 1;
        end
        
        plot(entropy,PSNR_score);
        figure = figure_num + 1;
    end
    show_result(conf_set);
end
%
