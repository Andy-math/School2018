function S = cal_s(Point)
global Bond_p Bond_num
all_pt = vertcat(Bond_p,Point);
temp = 1:Bond_num;
C = horzcat(transpose(temp),transpose(temp([2:end,1])));
DT = delaunayTriangulation(all_pt,C);
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

end