%inPoints = getPolygonGrid(lon,lat,ppa) returns points that are within a
%concave or convex polygon using the inpolygon function.

%lon and lat are columns representing the vertices of the polygon, as used
%in the Matlab function inpolygon

%ppa refers to the points per unit area you would like inside the polygon.
%Here unit area refers to a 1.0 X 1.0 square in the axes.

function [inPoints] = polygrid( object, ppa)
lon = []; lat = [];
for ii = 1:length(object.area)
    lon = [lon;object.area{ii}(:,2)];
    lat = [lat;object.area{ii}(:,1)];
end
N = sqrt(ppa);

% Find the bounding rectangle
lower_lon = min(lon);
higher_lon = max(lon);
lower_lat = min(lat);
higher_lat = max(lat);

% Create a grid of points within the bounding rectangle
inc_lon = 1/N;
inc_lat = 1/N;
interval_lon = lower_lon:inc_lon:higher_lon;
interval_lat = lower_lat:inc_lat:higher_lat;
[bigGridLon, bigGridLat] = meshgrid(interval_lon, interval_lat);

%Filter grid to get only points in polygons
inPoints = [];
for ii = 1:length(object.area)
    lon = object.area{ii}(:,2);
    lat = object.area{ii}(:,1);
    inallowed = inpolygon(bigGridLon(:), bigGridLat(:), lon, lat);
    for jj = 1:length(object.forbiddenarea)
        lonforbidden = object.forbiddenarea{jj}(:,2);
        latforbidden = object.forbiddenarea{jj}(:,1);
        inforbidden = inpolygon(bigGridLon(:), bigGridLat(:),...
            lonforbidden, latforbidden);
        for kk = 1:length(inallowed)
            if inforbidden(kk) == 1
                inallowed(kk) = 0;
            end
        end
    end
    inPoints = [inPoints;[bigGridLon(inallowed), bigGridLat(inallowed)]];
end

end