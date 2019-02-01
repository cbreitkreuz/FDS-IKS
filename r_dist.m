function dist = r_dist(x,y,xx,yy)
% Function r_dist computes distances from coordinates (x,y) = longitude,
% latitude to all coordinates in xx and yy


if x == 0 && y == 90
     warning('r_dist: x = 0 and y= 90!! For these values the function DOES NOT WORK!')
end

  pi180 = pi/180;
  
  earth_radius = 6378.137;


  lon1 = x*pi180;
  lon2 = xx*pi180;
  lat1 = y*pi180;
  lat2 = yy*pi180;

  dlon = lon2 - lon1;
  dlat = lat2 - lat1;
  a = (sin(dlat/2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon/2)).^2;
  
  temp =1-a;
  
  index = find(temp<0);
  temp(index) = 0;
  
  angles = 2 * atan2( sqrt(a), sqrt(temp) );
  dist = earth_radius * angles;

return