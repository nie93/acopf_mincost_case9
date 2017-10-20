%% ARCHIVE YALMIP
% yalmip solutions

% Parameters
c2 = mpc.gencost(:,5);
c1 = mpc.gencost(:,6);
c0 = mpc.gencost(:,7);
pd = mpc.bus(:,3);
qd = mpc.bus(:,4);


% Variables
e = sdpvar(numOfBus,1);
f = sdpvar(numOfBus,1);
pinj = sdpvar(numOfBus,1);
qinj = sdpvar(numOfBus,1);
pg = sdpvar(numOfGen,1);
qg = sdpvar(numOfGen,1);

% Constraints
constraint_yalmip = [];
for i = 1:numOfBuses
    constraint_yalmip = [constraint_yalmip, G(i,:)*(e(i)*e + f(i)*f) + B(i,:)*(f(i)*e - e(i)*f) == pinj(i)];
    constraint_yalmip = [constraint_yalmip, G(i,:)*(f(i)*e - e(i)*f) - B(i,:)*(e(i)*e + f(i)*f) == qinj(i)];
    constraint_yalmip = [constraint_yalmip, e(i)^2 + f(i)^2 <= vMagMax(i)^2];
    constraint_yalmip = [constraint_yalmip, e(i)^2 + f(i)^2 >= vMagMin(i)^2];
    if(any(busNumOfPV == i) || any(busNumOfSlack == i))
        constraint_yalmip = [constraint_yalmip, pg(i) - pd(i) == pinj(i)];
        constraint_yalmip = [constraint_yalmip, qg(i) - qd(i) == qinj(i)];
    end
    
    if(any(busNumOfPQ == i))
        constraint_yalmip = [constraint_yalmip, - pd(i) == pinj(i)];
        constraint_yalmip = [constraint_yalmip, - qd(i) == qinj(i)];
    end
    
end


obj_yalmip = pg'*diag(c2)*pg + c1'*pg + sum(c0);

assign(e,ones(numOfBus,1));
assign(f,zeros(numOfBus,1));
% assign(e,result_opf.bus(:,8).*cos(2*pi/360*result_opf.bus(:,9)));
% assign(f,result_opf.bus(:,8).*sin(2*pi/360*result_opf.bus(:,9)));
assign(pg,result_opf.gen(:,2));
ops = sdpsettings('usex0',1,'fmincon.MaxIter',10000);
optimize(constraint_yalmip,obj_yalmip,ops)
value(obj_yalmip)
[value(e), value(f)]


%%
% Variables
e = sdpvar(numOfBus,1);
f = sdpvar(numOfBus,1);
pg = sdpvar(numOfGen,1);
qg = sdpvar(numOfGen,1);

pinj([busNumOfSlack;busNumOfPV]) = pg - pd([busNumOfSlack;busNumOfPV]);
qinj([busNumOfSlack;busNumOfPV]) = qg - qd([busNumOfSlack;busNumOfPV]);
pinj(busNumOfPQ) = - pd(busNumOfPQ);
qinj(busNumOfPQ) = - qd(busNumOfPQ);

% Constraints
constraint_yalmip = [];
for i = 1:numOfBuses
    constraint_yalmip = [constraint_yalmip, G(i,:)*(e(i)*e + f(i)*f) + B(i,:)*(f(i)*e - e(i)*f) == pinj(i)];
    constraint_yalmip = [constraint_yalmip, G(i,:)*(f(i)*e - e(i)*f) - B(i,:)*(e(i)*e + f(i)*f) == qinj(i)];
    constraint_yalmip = [constraint_yalmip, e(i)^2 + f(i)^2 <= vMagMax(i)^2];
    constraint_yalmip = [constraint_yalmip, e(i)^2 + f(i)^2 >= vMagMin(i)^2];
end


obj_yalmip = pg'*diag(c2)*pg + c1'*pg + sum(c0);

% assign(e,ones(numOfBus,1));
% assign(f,zeros(numOfBus,1));
assign(e,result_opf.bus(:,8).*cos(2*pi/360*result_opf.bus(:,9)));
assign(f,result_opf.bus(:,8).*sin(2*pi/360*result_opf.bus(:,9)));
assign(pg,result_opf.gen(:,2));
ops = sdpsettings('usex0',1,'fmincon.MaxIter',10000);
optimize(constraint_yalmip,obj_yalmip,ops)
value(obj_yalmip)
[value(e), value(f)]


%%

v = sdpvar(numOfBus,1);
delta = sdpvar(numOfBus,1);
pg = sdpvar(numOfGen,1);
qg = sdpvar(numOfGen,1);

pinj([busNumOfSlack;busNumOfPV]) = pg - pd([busNumOfSlack;busNumOfPV]);
qinj([busNumOfSlack;busNumOfPV]) = qg - qd([busNumOfSlack;busNumOfPV]);
pinj(busNumOfPQ) = - pd(busNumOfPQ);
qinj(busNumOfPQ) = - qd(busNumOfPQ);

% Constraints
constraint_yalmip = [];
for i = 1:numOfBuses
    sumP = []; 
    sumQ = [];
    for j = 1:numOfBuses
        sumP = sumP + v(j)*( G(i,j)*cos(delta(i) - delta(j)) + B(i,j)*sin(delta(i) - delta(j)) );
        sumQ = sumQ + v(j)*( G(i,j)*cos(delta(i) - delta(j)) - B(i,j)*sin(delta(i) - delta(j)) );
    end
    
    constraint_yalmip = [constraint_yalmip, v(i).*sumP == pinj(i)];
    constraint_yalmip = [constraint_yalmip, v(i).*sumQ == qinj(i)];
    constraint_yalmip = [constraint_yalmip, v(i) <= vMagMax(i)];
    constraint_yalmip = [constraint_yalmip, v(i) >= vMagMin(i)];
end
