function points = PolygonReduction(pointList, epsilon, count)

dmax = 0;
index = 0;

resultList = 0;

l = length(pointList);

startPoint = pointList(1, :);

if(count == 0)
    [endPoint, distance] = FarthestPoint(startPoint, pointList(2:end, :));
else
    endPoint = pointList(end, :);
end

for i = 2:(l-1)
    d = Distance2Line(startPoint, endPoint,  pointList(i, :));
    if( d > dmax )
        index = i;
        dmax = d;
    end
end

count = count + 1;

if(dmax > epsilon)
    
    recResult1 = PolygonReduction(pointList(1:index, :), epsilon, count);   
    recResult2 = PolygonReduction(pointList(index:end, :), epsilon, count);
    
    resultList = [recResult1(1:(end-1), :); recResult2(:, :)];
else
    resultList = [pointList(1, :); pointList(end, :)];
end

points = resultList;

end