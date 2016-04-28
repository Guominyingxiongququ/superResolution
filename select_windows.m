%% f and g are two image sets 
%  this function is used to select the patches 
%  which has the worst performance based on PSNR
%  select the worst 20% patch
function window_list = select_windows(conf)
    for j = 3
        results_truth = cell(numel(conf.results), 1);
        results_test = cell(numel(conf.results), 1);
        scale = 1;
        for i = 1:numel(conf.results)
            results_truth{i} = conf.results{i}{1};
            results_test{i} = conf.results{i}{j};
        end
        [img_truth,~,~] = load_images(results_truth);
        [img_test,~,~] = load_images(results_test);
        conf.window = [99,99]
        window_truth = collect_patches(conf, img_truth, scale); 
        window_test = collect_patches(conf, img_test, scale);
        select_num = ceil(size(window_truth,3)/5);
        window_list = select_worst_patch(window_truth,window_test,select_num);
    end
end

function patch_list = select_worst_patch(patches_truth,patches_test,select_num) 
    score = zeros(1,size(patches_truth,3));
    for i = 1:size(patches_truth,3)
        E = patches_truth(:,:,i) - patches_test(:,:,i);
        N= numel(E);
        score(1,i) = 10*log10( N / sum(E(:).^2) );
    end
    [sortedValue_score , score_Ranked] = sort(score,'ascend'); 
    patch_list = zeros(size(patches_truth,1),size(patches_truth,2),select_num);
    for i = 1:select_num
        patch_list(:,:,i) = patches_truth(:,:,score_Ranked(i));
    end
end