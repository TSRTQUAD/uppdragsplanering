%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Missionplaner Version 1.0             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
startpoint = [15.5737 58.3948];
endpoint = [15.573705 58.3948];
load('object.mat')
normalizationfactor = inf;

% ----------------- Place nodes --------------------
% Create nodes within the polygon
inPoints = polygrid(object,4e8);

% --------------- Find costmatrix ------------------
% Create cost matrixes, 3 different fly patterns
tmpcmat = createcostmatrix(inPoints,object);

% The TSP-solver kan only handle integers
normalizationfactor = min(2^16/max([max(max(tmpcmat.lon)),max(max...
    (tmpcmat.lat)),max(max(tmpcmat.cross))]),normalizationfactor);

% --------------- Find trajectory ------------------
pattern = {'lon' 'lat' 'cross'}; finaltrajectory = [];
nrofnodes = length(inPoints);
bestsolution = inf;
for jj = 1:3
    cmat = round(normalizationfactor*tmpcmat.(pattern{jj}));
    trajectory.(pattern{jj}) = tspsolver(inPoints,cmat);
    
    % Calculate the total trajectory distance and sort out the shortest
    distance = 0;
    for kk = 2:length(trajectory.(pattern{jj}))
        distance = distance + lldistkm(trajectory...
            .(pattern{jj})(kk-1,:),trajectory.(pattern{jj})(kk,:));
    end
    currentsolution = distance;
    if currentsolution < bestsolution
        besttrajectory = trajectory.(pattern{jj});
        bestsolution = currentsolution;
    end
end
rawtrajectory = [startpoint; besttrajectory; endpoint];
% ------------ Interpolate trajectory --------------
% Interpolate using parametric splines
trajectory = interparc(2e3,rawtrajectory(:,1),rawtrajectory(:,2),'spline');

% --------------- Present results ------------------
% Calculate the total trajectory distance
distance = 0;
for ii = 2:length(trajectory)
    distance = distance + lldistkm(trajectory(ii-1,:),trajectory(ii,:));
end

% Plot and calculate the covered area in m^2
earthellipsoid = referenceSphere('earth','m');
earthellipsoidsurfacearea = areaquad(-90,-180,90,180,earthellipsoid);
figure('Name','Optimal trajectory for area coverage','Numbertitle','off')
clf; hold on; area = 0; forbiddenarea = 0;
for ii = 1:length(object.area)
    area = area + areaint(object.area{ii}(:,1),object.area{ii}(:,2))*...
        earthellipsoidsurfacearea;
    h1 = fill(object.area{ii}(:,2),object.area{ii}(:,1),[0.5,0.5,0.5]);
end
for ii = 1:length(object.forbiddenarea)
    forbiddenarea = forbiddenarea+areaint(object.forbiddenarea{ii}(:,1),...
        object.forbiddenarea{ii}(:,2))*earthellipsoidsurfacearea;
    fill(object.forbiddenarea{ii}(:,2),object.forbiddenarea{ii}(:,1),'w');
end
area = area - forbiddenarea;
plot(inPoints(:,1),inPoints(:,2), '.k');
h2 = plot(trajectory(:,1),trajectory(:,2),'r');
legend([h1 h2],{['Total area search: ' num2str(area) 'm^2'],...
    ['Total trajectory distance: ' num2str(distance*1000) 'm']})
hold off

