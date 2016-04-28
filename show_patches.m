function show_patches(patch_list, conf,row_size,scale,figure_num)
    [~,~,patch_num] = size(patch_list);
    row_num = ceil(patch_num/row_size);
    patch_row = cell(1,row_num);
    zero_patch = zeros(conf.window(1)*scale,conf.window(2)*scale);
    for i = 1:row_num
        for j = 1:row_size
            if (row_size*(i-1)+j)<=patch_num
                if j==1
                    cur_patch = patch_list(:,:,row_size*(i-1)+1);
                else
                    cur_patch = cat(2,cur_patch,patch_list(:,:,row_size*(i-1)+j));
                end
            else
                cur_patch = cat(2,cur_patch,zero_patch);
            end
        end
        patch_row{i} = cur_patch;
        if i==1
            complete_patch = patch_row{1};
        else
            complete_patch =cat(1, complete_patch,patch_row{i}); 
        end
    end
    figure(figure_num);
    imshow(complete_patch);
end