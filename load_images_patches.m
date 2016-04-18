function [imgs imgsCB imgsCR] = load_images_patches(paths,conf)

patch_size = conf.patch_size(1);
imgs = cell(0);
imgsCB = cell(0);
imgsCR = cell(0);
num = 0;
for i = 1:numel(paths)
    X = imread(paths{i});
        X = rgb2ycbcr(X);        
        [length,width,~] = size(X);
        x_patch_num = floor(length/patch_size);
        y_patch_num = floor(width/patch_size);
        if (x_patch_num ~= 0)&&(y_patch_num ~= 0)
            for j =1:x_patch_num
                for k = 1:y_patch_num
                    num = num + 1;
                    curImgs = cell(1);
                    %X = rgb2gray(X);      
                    if size(X, 3) == 3 % we extract our features from Y channel
                        offset_x = (j-1) * patch_size+1;
                        offset_y = (k-1) * patch_size+1;
                        curImgsCB = cell(1);
                        curImgsCR = cell(1);
                        curImgsCB{1} = im2single(X(offset_x:offset_x+patch_size-1,offset_y:offset_y+patch_size-1,2)); 
                        curImgsCR{1} = im2single(X(offset_x:offset_x+patch_size-1,offset_y:offset_y+patch_size-1,3));
                        X = X(offset_x:offset_x+patch_size-1,offset_y:offset_y+patch_size-1,1);
                    end
                    X = im2single(X); % to reduce memory usage
                    curImgs{1} = X;
                    imgs = [imgs curImgs];
                    imgsCB = [imgsCB curImgsCB];
                    imgsCR = [imgsCR curImgsCR];
                end
            end
        end
end
