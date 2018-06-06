clear
Bond_p = [-2 0;-2 2;-1 2;-1 1;1 1;1 2;2 2;2 0];
[Bond_num,~] = size(Bond_p);
pgon = polyshape(Bond_p);
figure(1)
subplot(131);plot(pgon)
Point_Num = 200;
rng shuffle
Point(:,1) = rand(Point_Num,1)*4-2;
Point(:,2) = rand(Point_Num,1)*2;
subplot(132);scatter(Point(:,1),Point(:,2))
TFin = isinterior(pgon,Point);
Point = Point(TFin,:);
subplot(133);scatter(Point(:,1),Point(:,2))
all_pt = vertcat(Bond_p,Point);
temp = 1:Bond_num;
global C
C = horzcat(transpose(temp),transpose(temp([2:end,1])));
DT = delaunayTriangulation(all_pt,C);
[K,~] = convexHull(DT);
figure(2)
axis([-3 3 -1 3])
IO = isInterior(DT);
triplot(DT(IO, :),DT.Points(:,1), DT.Points(:,2),'LineWidth', 1)
usf_cn = DT.ConnectivityList(IO,:);
[row,~] = size(usf_cn);
convx_p = reshape(DT.Points(usf_cn,:),[],3,2);
convx_eg(:,1) = sqrt(sum((convx_p(:,1,:)-convx_p(:,2,:)).^2,3));
convx_eg(:,2) = sqrt(sum((convx_p(:,1,:)-convx_p(:,3,:)).^2,3));
convx_eg(:,3) = sqrt(sum((convx_p(:,2,:)-convx_p(:,3,:)).^2,3));
convx_eg(:,4) = sum(convx_eg,2)/2;
Area =sqrt(convx_eg(:,4).*...
    (convx_eg(:,4)-convx_eg(:,1)).*...
    (convx_eg(:,4)-convx_eg(:,2)).*...
    (convx_eg(:,4)-convx_eg(:,3))) ; % Heron's formula
ave = sum(Area)/row;
S = sum((Area-ave).^2);