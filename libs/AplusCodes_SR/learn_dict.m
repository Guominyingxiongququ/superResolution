function [conf] = learn_dict(conf, hires, dictsize)
% Sample patches (from high-res. images) and extract features (from low-res.)
% for the Super Resolution algorithm training phase, using specified scale 
% factor between high-res. and low-res.

% Load training high-res. image set and resample it
temp_patches = hires;
hires = modcrop(hires, conf.scale); % crop a bit (to simplify scaling issues)
% Scale down images
lores = resize(hires, 1/conf.scale, conf.interpolate_kernel);

midres = resize(lores, conf.upsample_factor, conf.interpolate_kernel);
features = collect(conf, midres, conf.upsample_factor, conf.filters);
clear midres

interpolated = resize(lores, conf.scale, conf.interpolate_kernel);
clear lores
patches = cell(size(hires));
for i = 1:numel(patches) % Remove low frequencies
    patches{i} = hires{i} - interpolated{i};
end
% clear hires interpolated

patches = collect(conf, patches, conf.scale, {});

% feature size * number of features
% Set KSVD configuration
%ksvd_conf.iternum = 20; % TBD
ksvd_conf.iternum = 20; % TBD
ksvd_conf.memusage = 'normal'; % higher usage doesn't fit...
%ksvd_conf.dictsize = 5000; % TBD
ksvd_conf.dictsize = dictsize; % TBD
ksvd_conf.Tdata = 3; % maximal sparsity: TBD
ksvd_conf.samples = size(patches,2);

ksvd_window_conf = ksvd_conf;

% PCA dimensionality reduction
C = double(features * features');
[V, D] = eig(C);
D = diag(D); % perform PCA on features matrix 
D = cumsum(D) / sum(D);
k = find(D >= 1e-3, 1); % ignore 0.1% energy
conf.V_pca = V(:, k:end); % choose the largest eigenvectors' projection
features_pca = conf.V_pca' * features;
conf.patch_num = 1000;
%k means
numClusters = ceil(size(features,2)/conf.cluster_size);

% use kmeans

if conf.kmeans_window == 1
    [centers1,index] = vl_kmeans(features_pca, conf.word_num, 'Algorithm', 'ANN');
    [update_centers1 , center_index1] = find_closest_centers(centers1 , index, features_pca);  
    windows_hist = zeros(conf.word_num,numel(hires));
    window_patch_size = size(patches,2)/numel(hires); %record how many patches each window has
    for i = 1:numel(hires)
        for j = 1:window_patch_size
            word = index(1, (i-1)*window_patch_size + j);
            windows_hist(word,i)=windows_hist(word,i)+1;
        end
    end
    cluster_num = ceil(numel(hires)/4);
    [centers, index] = vl_kmeans(windows_hist,cluster_num, 'Algorithm', 'ANN', 'MaxNumComparisons', 1000);
    [update_centers , center_index] = find_closest_centers(centers , index, windows_hist);
    update_center_num = size(update_centers,2);
    update_feature_size =  update_center_num * window_patch_size;
    update_features = zeros(size(features_pca,1),update_feature_size);
    for i =1:update_center_num
        offset = window_patch_size;
        update_begin= (i-1)*window_patch_size+1;
        begin= (center_index(1,i)-1)*window_patch_size+1;
        update_features(:,update_begin:update_begin+offset-1) = features_pca(:,begin:begin+offset-1);
        update_patches(:,update_begin:update_begin+offset-1) = patches(:,begin:begin+offset-1);
    end
    features_pca = update_features;
    patches = update_patches;
end

% <<<<<<< HEAD
% 
% =======
% %% Set KSVD configuration
% %ksvd_conf.iternum = 20; % TBD
% ksvd_conf.iternum = 20; % TBD
% ksvd_conf.memusage = 'normal'; % higher usage doesn't fit...
% ksvd_conf.dictsize = dictsize; % TBD
% ksvd_conf.Tdata = 3; % maximal sparsity: TBD
% ksvd_conf.samples = size(patches,2);
% conf.ksvd_conf = ksvd_conf;
% 
% %%
% >>>>>>> origin/master
if conf.kmeans == 1 
    [centers,index] = vl_kmeans(features_pca, 1000, 'Algorithm', 'ANN', 'Initialization','RANDSEL');
    u_index = unique(index);
    MAX = 100000;
    newCenters = zeros(size(centers));

    closestPatchC = zeros(1,size(centers,2));
    center_distance = ones(1,size(centers,2))* MAX;
    newPatches = zeros(size(patches,1),size(centers,2));
    
    for i = 1:size(index,2)
        c = index(1,i);  %find the center corresponded to the patch i
        dis = getDistance(centers(:,c),features_pca(:,i));
        if dis < center_distance(1,c)
            center_distance(1,c) = dis;
            newCenters(:,c) = features_pca(:,i);
            closestPatchC(1,c) = i;
            newPatches(:,c) = patches(:,i);
        end
    end
else
    if conf.ksvd == 1
        ksvd_conf.data = double(features_pca);
        fprintf('Training [%d x %d] dictionary on %d vectors using K-SVD\n', ...
        size(ksvd_conf.data, 1), ksvd_conf.dictsize, size(ksvd_conf.data, 2))
        [conf.dict_lores, gamma] = ksvd(ksvd_conf);
    else
        select_num = conf.patch_num;
        total_num = size(features_pca,2);
        select_patch_list = sort(randperm(total_num,select_num));
    end
end
 
clear C D V

if conf.kmeans == 1
    conf.dict_lores = newCenters(:,u_index);
    conf.dict_hires = newPatches(:,u_index);
else
    if conf.ksvd == 1
        conf.dict_hires = (patches * (full(gamma))') * inv(full(gamma * gamma'));
    else
        conf.dict_lores = features_pca(:,select_patch_list);
        conf.dict_hires = patches(:,select_patch_list);
    end
end
clear features_pca


fprintf('Computing high-res. dictionary from low-res. dictionary\n');
% patches = double(patches(:,index)); % Since it is saved in single-precision.
% conf.dict_hires = double(dict_hires); 