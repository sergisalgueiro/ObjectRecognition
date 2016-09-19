function dis = Distance2Line(P1, P2, p)
%https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points
xdif = (P2(1) - P1(1));
ydif = (P2(2) - P1(2));

gap = sqrt(ydif^2 + xdif^2);

num = abs(ydif*p(1) - xdif*p(2) + P2(1)*P1(2) - P2(2)*P1(1));

dis = num/gap;

end