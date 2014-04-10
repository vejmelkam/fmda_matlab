

function gcd = great_circle_distance(lat1,lon1,lat2,lon2)

    rlat1 = pi * lat1 / 180.0;
    rlat2 = pi * lat2 / 180.0;
    rlon1 = pi * lon1 / 180.0;
    rlon2 = pi * lon2 / 180.0;
    
    a = sin(0.5*(rlat1 - rlat2))^2 + cos(rlat1)*cos(rlat2)*sin(0.5*(rlon1 - rlon2))^2;
    c = 2 * atan2(a^0.5, (1-a)^0.5);
    gcd = 6371.0 * c;
    
end