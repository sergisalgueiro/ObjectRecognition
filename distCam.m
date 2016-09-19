function [distance, front, side, angle] = distCam(tbBox, scale, height, dImg)

    focal_mm=3.6;
    sensor_width_mm=3.76;
    image_width_in_pixels=1024;
    focalLength = (focal_mm / sensor_width_mm) * image_width_in_pixels;
%     focalLength = 455;
    
    objectSizePix =  tbBox(4) * scale^-1;
    objectSizeReal = height;
    distance = focalLength * objectSizeReal / objectSizePix;

%     horizontalPix = (tbBox(1) + tbBox(3)/2) * scale^-1;
%     imgWidthPix = image_width_in_pixels;
    angImg = -((tbBox(1) + tbBox(3)/2)-(size(dImg,2)/2))*(53.5/2)/(size(dImg,2)/2);
    %negative because of robot Odometry plane
    front = distance*cos(deg2rad(angImg));
    side = distance*sin(deg2rad(angImg));
    angle=angImg;
    
%     if angle >0
%         angle = pi/2-angle;
%     elseif angle < 0
%         angle = -pi/2 - angle;
%     end
end