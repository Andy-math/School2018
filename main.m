clear
Point_Num = 200;
z_expd = 9;
Point = zeros(Point_Num,3);
rng('shuffle')
Point(:,1) = rand(Point_Num,1)*2*pi;
Point(:,2) = rand(Point_Num,1);
Point(1:Point_Num/2,3) = sqrt(z_expd-z_expd*Point(1:Point_Num/2,2).^2);
Point(Point_Num/2+1:Point_Num,3) = -sqrt(3-3*Point(Point_Num/2+1:Point_Num,2).^2);
[Point(:,1),Point(:,2)] = pol2cart(Point(:,1),Point(:,2));

options.Algorithm = 'active-set';
options.Display = 'iter-detailed';
options.MaxFunctionEvaluations = inf;
options.MaxIterations = inf;
options.UseParallel = true;
Point = fmincon(@cal_s,Point,[],[],[],[],[],[],@mycon,options);

function [c,ceq] = mycon(Point)
    c = [];
    ceq = Point(:,1).^2+Point(:,2).^2+Point(:,3).^2/9-1;
end

