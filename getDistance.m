function dis = getDistance(vector1, vector2)
    length = size(vector1,1);
    dis = 0;
    for i = 1:length        
        dis=dis+(vector1(i,1) - vector2(i,1))^2;
    end
    dis = sqrt(dis);
end