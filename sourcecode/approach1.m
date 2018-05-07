%read data in
dt = importdata('C:\Users\weilan\Documents\513\final_project_data\final_project_data\final_project_point_cloud.fuse', ' ');

a=dt(:,1);
b=dt(:,2);
c=dt(:,3);
d=dt(:,4);
dt3 = [a b c];
%lla to ecef
dtECEF = lla2ecef(dt3,'WGS84');

lat0 = 45.90414414;
lon0 = 11.02845385;
%ecef to ned
[uNorth,vEast,wDown] = ecef2nedv(dtECEF(:,1),dtECEF(:,2),dtECEF(:,3),lat0,lon0);

Qs = 0.362114;
Qx = 0.374050;
Qy = 0.592222;
Qz = 0.615007;
Rq = [[Qs*Qs+Qx*Qx-Qy*Qy-Qz*Qz 2*Qx*Qy-2*Qs*Qz 2*Qx*Qz+2*Qs*Qy];[2*Qx*Qy+2*Qs*Qz Qs*Qs-Qx*Qx+Qy*Qy-Qz*Qz 2*Qy*Qz-2*Qs*Qx];[2*Qx*Qz-2*Qs*Qy 2*Qz*Qy+2*Qs*Qx Qs*Qs-Qx*Qx-Qy*Qy+Qz*Qz]];
%end to camera coordinates
dtCC = [uNorth vEast wDown]*Rq;

pcshow(dtCC, d, 'MarkerSize', 1000);
view([180,90]);

dtFull = [dtCC d];
%select the high altitute section on image
edgeLimi = 0.25*(max(dtFull(:,2))-min(dtFull(:,2)))+min(dtFull(:,2));
dtHigh = dtFull(dtFull(:,2)<edgeLimi, :);

dtRest = setdiff(dtFull,dtHigh,'rows');
dtHighHL = [[dtHigh(:, [1,2,3]) 255*ones(length(dtHigh),1)];dtRest];

pcshow(dtHighHL(:, [1,2,3]), dtHighHL(:, 4), 'MarkerSize', 1000);
view([180,90]);
%select the lamp pole position
maxH = max(dtHigh(:,1))+2;
minH = min(dtHigh(:,1))-2;
dtKeep = dtFull(dtFull(:,1)<=maxH,:);
dtKeep = dtKeep(dtKeep(:,1)>=minH,:);
paraG1 = 0.48; % can be tuned 
paraG2 = 0.54; % can be tuned 
ground1 = paraG1*(max(dtFull(:,2))-min(dtFull(:,2)))+min(dtFull(:,2));
ground2 = paraG2*(max(dtFull(:,2))-min(dtFull(:,2)))+min(dtFull(:,2));
edge1 = dtKeep(dtKeep(:,1)<=maxH,:);
edge1 = edge1(edge1(:,1)>=maxH-5,:);
edge1 = edge1(edge1(:,2)<ground1,:);
edge2 = dtKeep(dtKeep(:,1)>=minH,:);
edge2 = edge2(edge2(:,1)<=minH+5,:);
edge2 = edge2(edge2(:,2)<ground2,:);
dtEdge = [edge1;edge2];

dtRestNew = setdiff(dtRest,dtEdge,'rows');
%final poles highlighted image
dtPoleHL = [[dtEdge(:,[1,2,3]) 255*ones(length(dtEdge),1)];[dtHigh(:,[1,2,3]) 255*ones(length(dtHigh),1)];dtRestNew];

pcshow(dtPoleHL(:, [1,2,3]), dtPoleHL(:, 4), 'MarkerSize', 1000);
view([180,90]);

%save data for approach II
dlmwrite('matlabData.csv', dtPoleHL(:, [1,2,3]), 'precision', '%10.16f'); 