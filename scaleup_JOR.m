% image scale-up via JOR
% by dengxin dai
% Aug 27 2014

function [imgs, midres] = scaleup_JOR(conf, JOR, imgs)

% Super-Resolution Iteration
    fprintf('Scale-Up via JOR ... ');
    midres = resize(imgs, conf.upsample_factor, conf.interpolate_kernel);    
    
    for i = 1:numel(midres)
        
        %% extract features
        im_interpolated = double(midres{i});
        
        % gaurantee  to use the same pca as training preference forest    
        features = collect(conf, {im_interpolated}, conf.scale, conf.filters);
%         features = double(features);          
        features = single(JOR.V_pca')*features;

        % normalize features
        l2 = sum(features.^2).^0.5+eps;
        l2n = repmat(l2,size(features,1),1);   
        features_knn = features./l2n;       
    
        %% vl_feat 
        nn = JOR.knn_k; img_num = size(features, 2); 
        [nnind, distance] = vl_kdtreequery(JOR.kdtree, JOR.lowpatches, single(features_knn), 'NumNeighbors', nn,'MaxComparisons', JOR.kdtree_hitnum) ;
       
        
        % error table
        error_temp = reshape(JOR.error_table(nnind,:), [nn img_num JOR.reg_num]);

        % weighted 
        weights = nn:-1:1; weights = repmat((weights(:)), [1, img_num JOR.reg_num]); 
        nn_error = sum(error_temp .* weights, 1);
     
%         nn_error = sum(error_temp);
        
        nn_error = squeeze(nn_error); 

        
        [~, plabel_nn] =  min(nn_error, [], 2);   

        for l = 1:size(features,2)            
            patches(:, l) = JOR.PPs{plabel_nn(l)} * features(:,l);
%                     recon_error(l,1) = nn_error(l,plabel_nn(l));
        end
           
 
       % Add low frequencies to each reconstructed patch        
        patches = patches + collect(conf, {im_interpolated}, conf.scale, {});         
      
        % Combine all patches into one image
        img_size = size(imgs{i}) * conf.scale;
        grid = sampling_grid(img_size, ...
            conf.window, conf.overlap, conf.border, conf.scale);
        result = overlap_add(patches, img_size, grid);                   

        imgs{i} = result; % for the next iteration
        fprintf('.');
    end
fprintf('JOR finished \n');




   

