

function [lati,loni,dist] = find_closest_point(lat,lon,glat,glon)

    dists = great_circle_distance(lat,lon,glat,glon);
    
    [dist,i] = min(dists(:));
    [loni,lati] = ind2sub(size(dists),i);
    