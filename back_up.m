clear
Point_Num = 100;
z_expd = 9;
Point = zeros(Point_Num,3);
rng('shuffle')
Point(:,1) = rand(Point_Num,1)*2*pi;
Point(:,2) = rand(Point_Num,1);
Point(1:Point_Num/2,3) = sqrt(z_expd-z_expd*Point(1:Point_Num/2,2).^2);
Point(Point_Num/2+1:Point_Num,3) = -sqrt(3-3*Point(Point_Num/2+1:Point_Num,2).^2);
[Point(:,1),Point(:,2)] = pol2cart(Point(:,1),Point(:,2));
DT = delaunayTriangulation(Point);
[K,v] = convexHull(DT);
trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3))

[row,~] = size(K);
convx_p = reshape(DT.Points(K,:),[],3,3);
convx_eg(:,1) = sqrt(sum((convx_p(:,1,:)-convx_p(:,2,:)).^2,3));
convx_eg(:,2) = sqrt(sum((convx_p(:,1,:)-convx_p(:,3,:)).^2,3));
convx_eg(:,3) = sqrt(sum((convx_p(:,2,:)-convx_p(:,3,:)).^2,3));
convx_eg(:,4) = sum(convx_eg,2)/2;
Area =sqrt(convx_eg(:,4).*...
    (convx_eg(:,4)-convx_eg(:,1)).*...
    (convx_eg(:,4)-convx_eg(:,2)).*...
    (convx_eg(:,4)-convx_eg(:,3))) ; % Heron's formula
% p=(a+b+c)/2
% S =sqrt[p(p-a)(p-b)(p-c)]
ave = sum(Area)/row;
S = sum((Area-ave).^2);
