function [patches] = collect_patches(conf, imgs, scale)

num_of_imgs = numel(imgs);
patch_cell = cell(num_of_imgs, 1); 
num_of_patches = 0;
for i = 1:num_of_imgs
    sz = size(imgs{i});
    grid = sampling_grid(sz, conf.window, conf.overlap, conf.border, scale);
    cur_patches = imgs{i}(grid);
    num_of_patches = num_of_patches + size(cur_patches,3);
    patch_cell{i} = cur_patches;
    patch_size = size(cur_patches,1);
end
clear imgs;
patches = zeros([patch_size patch_size num_of_patches], 'single');
offset = 0;
for i = 1:num_of_imgs
    F = patch_cell{i};
    N = size(F, 3); 
    patches(:, :, (1:N) + offset) = F;
    offset = offset + N;
end
end
