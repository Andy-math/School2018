clear
global Bond_p Bond_num
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

options.Algorithm = 'active-set';
options.Display = 'iter-detailed';
options.MaxFunctionEvaluations = inf;
options.MaxIterations = inf;
options.UseParallel = true;
Point = fmincon(@cal_s,Point,[],[],[],[],[],[],[],options);

function [c,ceq] = mycon(Point)
    c = [];
    ceq = [];
end