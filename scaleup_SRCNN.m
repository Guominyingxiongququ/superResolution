function [imgs, midres] = scaleup_SRCNN(conf, imgs)

% Super-Resolution Iteration
    fprintf('Scale-Up SRCNN ...');
    model = ['x' num2str(conf.scale) '.mat'];
    
    midres = resize(imgs, conf.upsample_factor, conf.interpolate_kernel);    
    
    for i = 1:numel(midres)
        
        
        % Combine all patches into one image
        img_size = size(imgs{i}) * conf.scale;
        im_h = SRCNN(model, double(midres{i}));

        imgs{i} = im_h; % for the next iteration
        fprintf('.');
    end
fprintf('\n');
