%% f and g are two image sets 
%  this function is used to select the patches 
%  which has the worst performance based on PSNR

function patch_list = select_patches(conf)
    for j = 2:numel(conf.results{1})
        results_truth = cell(numel(conf.results), 1);
        results_test = cell(numel(conf.results), 1);
        scale = 3;
        for i = 1:numel(conf.results)
            results_truth{i} = conf.results{i}{1};
            results_test{i} = conf.results{i}{j};
            [img_truth,~,~] = load_images(results_truth);
            [img_test,~,~] = load_images(results_test);
            patches_truth = collect_patches(conf, img_truth, scale); 
            patches_test = collect_patches(conf, imgs_test, scale);
            select_num = ceil(size(patches_truth,2)/5);
            patch_list = select_worst_patch(patches_truth,patches_test,select_num);
%             F = im2double(imread(f)); % original
%             G = im2double(imread(g)); % distorted
%             E = F - G; % error signal
%             N = numel(E); % Assume the original signal is at peak (|F|=1)
%             res = 10*log10( N / sum(E(:).^2) );
%             scores(i,j) = score;
        end
        
        
    end

end

function patch_list = select_worst_patch(patches_truth,patches_test,select_num) 
    score = zeros(1,size(patches_truth,2));
    for i = 1:size(patches_truth,2)
        E = patches_truth(:,i) - patches_test(:,i);
        N= numel(E);
        score(1,i) = 10*log10( N / sum(E(:).^2) );
    end
    [sortedValue_score , score_Ranked] = sort(score,'descend'); %%可能会有问题
    patch_list = zeros(size(patches_truth,1),select_num);
    for i = 1:select_num
        patch_list(:,i) = patches_truth(:,score_Ranked);
    end
end
% 
% %% collect file names from result 
% function result  = collect_fileName(c)
% 
% end
% 
% function res = calc_PeakSNR(f, g)
% F = im2double(imread(f)); % original
% G = im2double(imread(g)); % distorted
% E = F - G; % error signal
% N = numel(E); % Assume the original signal is at peak (|F|=1)
% res = 10*log10( N / sum(E(:).^2) );
% 
