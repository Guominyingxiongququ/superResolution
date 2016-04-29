function avg_entropy = get_entropy(window_list)
    window_num = size(window_list,3);
    sum_entropy = 0;
    for i = 1:window_num
        sum_entropy = sum_entropy + entropy(window_list(:,:,i));
    end
    avg_entropy = avg_entropy/window_num;
end