function [point, dis] = FarthestPoint(startPoint, points)

x = startPoint(1);
y = startPoint(2);

d = [];
dmax = 0;
index = 0;

for i = 1:length(points)
    
    if(points(i, 1) == x && points(i, 2) == y)
        continue;
    end
    
    d = sqrt((points(i, 1) - x)^2 + (points(i, 2) - y)^2);
    
    if( d > dmax )
        dmax = d;
        index = i;
    end
    
    
end

dis = dmax;
point = points(index, :);

end