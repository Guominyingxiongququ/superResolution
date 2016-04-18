% compute the reconstruction error of all patches for all projection
% functions, PPs is a set of projection funcions 

function [dict_member, recon_errors] = ProjectionErrorComput(lowpatches, highpatches, PPs)

    
    pps_num = numel(PPs);
        %membership 
       
% % % %         for j=1:pps_num   
% % % %             patch_recon = PPs{j}*lowpatches;
% % % %             recon_errors(:,j) = (mean( (highpatches - patch_recon).^2,  1)).^.5;      
% % % %         end
% % % %         
        
        
    range = 200000;   % to avoid out-of-memory 
    img_num = size(lowpatches, 2);
    
    recon_errors = zeros(img_num, pps_num);
    
    if img_num>range 
       
        iter_num = floor(img_num/range);
        for ii=1:iter_num
            q_ind = 1+(ii-1)*range : ii*range; 
            for j=1:pps_num   
                patch_recon = PPs{j}*lowpatches(:, q_ind);
                recon_errors(q_ind,j) = (mean( (highpatches(:, q_ind) - patch_recon).^2,  1)).^.5;      
            end
%             [ii iter_num]
        end

        if ii*range<img_num
            q_ind = ii*range+1:img_num; 
            for j=1:pps_num   
                patch_recon = PPs{j}*lowpatches(:, q_ind);
                recon_errors(q_ind,j) = (mean( (highpatches(:, q_ind) - patch_recon).^2,  1)).^.5;      
            end
        end 
     
        
    else 
        
        for j=1:pps_num   
            patch_recon = PPs{j}*lowpatches;
            recon_errors(:,j) = (mean( (highpatches - patch_recon).^2,  1)).^.5;      
        end
    end 
     
   [~, dict_member] = min(recon_errors, [], 2);
    
end
