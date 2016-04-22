function [update_centers , center_index] = find_closest_centers(centers , index, features)
    u_index = unique(index);
    
    MAX = 100000;
    update_centers = zeros(size(centers));
    center_index = zeros(1,size(centers,2));
    center_distance = ones(1,size(centers,2));
    center_distance = center_distance * MAX;
    
    for i = 1:size(index,2)
        c = index(1,i);  %find the center corresponded to the patch i
        dis = getDistance(centers(:,c),features(:,i));
        if dis < center_distance(1,c)
            center_distance(1,c) = dis;
            update_centers(:,c) = features(:,i);
            center_index(1,c) = i;
        end
    end
    
    update_centers = update_centers(:,u_index);
    center_index = center_index(:,u_index);
end