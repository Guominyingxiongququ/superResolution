% return mined dict_size projection functions PPS, membership of all samples dict_member, and the reconstruction errors
% of all samples by all functions 
function [PPs, dict_member, reerror_i] = Learn_Projections_EM(lowpatches, highpatches, dict_size, em_inum)        

       img_num = size(lowpatches, 2);
       
%% kmeans or random initialization        
%        lowpatches_s = lowpatches; 
%        if img_num > max_num   
%            lowpatches_s = lowpatches(:, 1:round(img_num/max_num):end); 
%        end 
%        [centervalues, kcenters] = kmeans(lowpatches_s, dict_size); 
%        [~, kcenters] = knn_classify(lowpatches', centervalues', 1:dict_size, 1); 
       
       kcenters =  randi(dict_size, [img_num, 1]); 
       
%         dict_size = conf.splits(split_id);
        num_limit = 4000;                     % the model is learned from this many patches (subsampled from all patches in the same class)
        for i = 1:dict_size
            idx = find(kcenters==i);
            if numel(idx)>num_limit        idx=idx(1:round(numel(idx)/num_limit):end); end
%             if numel(idx)>num_limit        [xxx, idx] = knn(centervalues(:,i), lowpatches, num_limit); end
            Lo = lowpatches(:, idx);       
%             PPs{i, 1} = highpatches(:,idx)*inv(Lo'*Lo+0.1*eye(size(Lo,2)))*Lo'; 
            PPs{i, 1} = highpatches(:,idx) * ( (Lo'*Lo+0.1*eye(size(Lo,2))) \ Lo' ); 
            
%             Hi*inv(Lo'*Lo+llambda*eye(size(Lo,2)))*Lo';
        end

        %membership 
        p_nums = size(lowpatches, 2);  
        [dict_member, reerror_i]= ProjectionErrorComput(lowpatches, highpatches, PPs);
    
        % multiple lambdas can be tried to get even better results
        lambdas = [0.15]; %[0.01 0.015 0.02 0.02 0.05 0.07 0.1 0.2 0.2 0.5 0.5 1 1.5 2 3 5 10 0.01 0.015 0.02 0.02 0.05 0.07 0.1 0.2 0.2 0.5 0.5 1 1.5 2 3 5 10]; 
            
        %iteration 
        if nargin< 4 em_inum = 20; end
        
        for em_i=1:em_inum

            % e-step: update PPs
            for i = 1:dict_size

                % update lambda 
                idx_all = find(dict_member==i);
                for lam_id=1:length(lambdas)           % try different lambdas  
                    
                    idx = find(dict_member==i);
                    if numel(idx)>num_limit idx = randperm(length(idx)); idx = idx_all(idx(1:num_limit));   end % idx=idx(1:round(numel(idx)/num_limit):end);    end %idx=idx(1:round(numel(idx)/num_limit):end); end
                    Lo = lowpatches(:, idx);     
                    
                    lambda = lambdas(lam_id); 
%                     PPs_i{lam_id} = highpatches(:,idx)*inv(Lo'*Lo+lambda*eye(size(Lo,2)))*Lo';      
                    PPs_i{lam_id} = highpatches(:,idx) * ( (Lo'*Lo+lambda*eye(size(Lo,2))) \ Lo');    

                end
                
                 [~, reerror_iidd]= ProjectionErrorComput(lowpatches(:, idx_all), highpatches(:,idx_all), PPs_i);
                  recon_errors_i = mean(reerror_iidd, 1);
                
                [~, min_lam_id] = min(recon_errors_i, [], 2); 
                PPs{i} = PPs_i{min_lam_id}; 
            end
            
           

            % m_step: update membership 
            [dict_member, reerror_i]= ProjectionErrorComput(lowpatches, highpatches, PPs);
            
            all_recon_error(em_i) =  mean(min(reerror_i, [], 2));
            if em_i>1
                fprintf('em the %d th iteration, with error reduction %f \n', em_i, all_recon_error(em_i-1)-all_recon_error(em_i)); 
            else 
                fprintf('em the %d th iteration, with error  %f \n', em_i, all_recon_error(em_i)); 
            end
            
            
        end

        
end