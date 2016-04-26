function [patches] = collect_patches(conf, imgs, scale)

num_of_imgs = numel(imgs);
patch_cell = cell(num_of_imgs, 1); 
for i = 1:num_of_imgs
    sz = size(imgs{i});
    F = sampling_grid(sz, conf.window, conf.overlap, conf.border, scale);
    num_of_patches = num_of_patches + size(F,2);
    patch_cell{i} = F;
    patch_size = size(F,1);
end
clear imgs;
patches = zeros([patch_size num_of_patches], 'single');
offset = 0;
for i = 1:num_of_imgs
    F = patch_cell{i};
    N = size(F, 2); % number of features in current cell
    patches(:, (1:N) + offset) = F;
    offset = offset + N;
end
end
