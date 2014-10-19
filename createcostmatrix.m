function tmpcmat = createcostmatrix(inPoints,object)

% Fixed variables and punishes
patternpunish = 0.1;
saftyfactor = 1.1;

% Create a distance/cost and temporary matrix
nPoints = size(inPoints,1);
a = meshgrid(1:nPoints);
dmat = reshape(sqrt(sum((inPoints(a,:)-inPoints(a',:)).^2,2)),nPoints,nPoints);
entitystep = min(dmat(find(dmat>0)))*saftyfactor;
boundarymat = zeros(size(dmat));
dmatlon = dmat;
dmatlat = dmat;
dmatcross = dmat;

% Extract boundaries
allboundaries = [object.area object.forbiddenarea];
boundarylines = [];
for ii = 1:length(allboundaries);
    for kk = 2:length(allboundaries{ii})
        boundarylines = [boundarylines;[allboundaries{ii}(kk-1,2)...
            allboundaries{ii}(kk-1,1) allboundaries{ii}(kk,2)...
            allboundaries{ii}(kk,1)]];
    end
end
nrofboundarylines = length(boundarylines);

% The additive punish for crossing a boundry is ~ 100 times the lowest distance
boundarypunish = 1e2*entitystep;

for ii = 2:nPoints
    for jj = (ii+1):nPoints
        % --------- Determine flypattern ----------
        % Cheep to fly longitudal
        if (inPoints(ii,1) == inPoints(jj,1))...
                && (dmat(ii,jj) <= entitystep)
            dmatlon(ii,jj) = dmat(ii,jj)*patternpunish;
            dmatlon(jj,ii) = dmatlon(ii,jj);
        end
        % Cheep to fly lattitudal
        if (inPoints(ii,2) == inPoints(jj,2))...
                && (dmat(ii,jj) <= entitystep)
            dmatlat(ii,jj) = dmat(ii,jj)*patternpunish;
            dmatlat(jj,ii) = dmatlat(ii,jj);
        end
        % Cheep to fly crosswise
        if (inPoints(ii,1) ~= inPoints(jj,1))...
                && (inPoints(ii,2) ~= inPoints(jj,2))...
                && (dmat(ii,jj) <= entitystep*1.5)...
                && (dmat(ii,jj) > entitystep)
            dmatcross(ii,jj) = dmat(ii,jj)*patternpunish;
            dmatcross(jj,ii) = dmatcross(ii,jj);
        end
        
        % ---------- Check boundarylines ----------
        if (dmat(ii,jj) <= 10*entitystep)
            linebetweenpoints = [inPoints(ii,1) inPoints(ii,2)...
                inPoints(jj,1) inPoints(jj,2)];
            intersections = nan*ones(1,2*nrofboundarylines);
            for kk = 1:nrofboundarylines
                [x,y] = lineintersect(linebetweenpoints,boundarylines(kk,:));
                intersections((2*kk-1):(2*kk)) = [x y];
            end
            status = any(not(isnan(intersections)));
            if status
                boundarymat(ii,jj) = boundarypunish;
                boundarymat(jj,ii) = boundarypunish;
            end
        end
    end
end

tmpcmat = struct('lon',dmatlon+boundarymat,'lat',dmatlat+boundarymat,...
    'cross',dmatcross+boundarymat);

end


%-----------------------------------------------------------------------
function [x,y] = lineintersect(l1,l2)
% This function finds where two lines (2D) intersect
% Each line must be defined in a vector by two points
% [P1X P1Y P2X P2Y], you must provide two arguments
warning('off')

% Default values for x and y, in case of error these are the outputs
x=nan; y=nan;
Pint = intersect(l1,l2);

% When the lines are paralel there's x or y will be Inf
if (any(Pint==Inf))
    Pint=[nan nan];
end

% Check if the value NaN were due to vetrical or horizontal lines, rotate
% the lines 45 degrees and make the intersectiontest again. If
% intersection found those values are used since they don't serve any
% nurmerical purpose.
if isnan(Pint)
    l1rot=[1/sqrt(2) -1/sqrt(2);1/sqrt(2) 1/sqrt(2)]*[l1(1) l1(3);l1(2) l1(4)];
    l1=[l1rot(1,1) l1rot(2,1) l1rot(1,2) l1rot(2,2)];
    l2rot=[1/sqrt(2) -1/sqrt(2);1/sqrt(2) 1/sqrt(2)]*[l2(1) l2(3);l2(2) l2(4)];
    l2=[l2rot(1,1) l2rot(2,1) l2rot(1,2) l2rot(2,2)];
    Pint = intersect(l1,l2);
end

% Put the solution inside x and y
x=Pint(2);y=Pint(1);

% Find maximum and minimum values for the final test
l1minX=min([l1(1) l1(3)]);
l2minX=min([l2(1) l2(3)]);
l1minY=min([l1(2) l1(4)]);
l2minY=min([l2(2) l2(4)]);

l1maxX=max([l1(1) l1(3)]);
l2maxX=max([l2(1) l2(3)]);
l1maxY=max([l1(2) l1(4)]);
l2maxY=max([l2(2) l2(4)]);

% Test if the intersection is a point from the two lines because
% all the performed calculations where for infinite lines
if ((x<l1minX) || (x>l1maxX) || (y<l1minY) || (y>l1maxY) ||...
        (x<l2minX) || (x>l2maxX) || (y<l2minY) || (y>l2maxY) )
    x=nan; y=nan;
end

end

%-----------------------------------------------------------------------
function Pint = intersect(l1,l2)
% Find coordinates of intersection
ml1=(l1(4)-l1(2))/(l1(3)-l1(1));
ml2=(l2(4)-l2(2))/(l2(3)-l2(1));
bl1=l1(2)-ml1*l1(1);
bl2=l2(2)-ml2*l2(1);
b=[bl1 bl2]';
a=[1 -ml1; 1 -ml2];
Pint=a\b;
end
