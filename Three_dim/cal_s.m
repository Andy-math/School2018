function S = cal_s(Point)

DT = delaunayTriangulation(Point);
[K,~] = convexHull(DT);
trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3))
drawnow

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
end