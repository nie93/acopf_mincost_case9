clear; clc;

%% Get data from MatPower
mpc = loadcase('case9');
mpc = loadcase('case14');
% mpc = loadcase('case24_ieee_rts');
result_pf = runpf(mpc);
result_opf = runopf(mpc);

%% MatPower OPF Results
clc;
cg_pf = get_cost(mpc,result_pf.gen(:,2));
cost_total_pf = sum(cg_pf);

cg_opf = get_cost(mpc,result_opf.gen(:,2));
cost_total_opf = sum(cg_opf);

%% ACOPF (Minimum Cost)
% % ================================================ % %
% % obj. --- sum(cg)
% % s.t. --- ac power flow equations
% % ================================================ % %
numOfBus = size(mpc.bus,1);   
numOfGen = size(mpc.gen,1);   
Ybus = makeYbus(mpc);
G = real(Ybus);
B = imag(Ybus);
busNumOfSystem = mpc.bus(:,1);
busNumOfSlack = find(mpc.bus(:,2) == 3);
busNumOfPV = find(mpc.bus(:,2) == 2);
busNumOfPQ = find(mpc.bus(:,2) == 1);
numOfBuses = length(busNumOfSystem);
numOfPVBuses = length(busNumOfPV);
numOfPQBuses = length(busNumOfPQ);
voltMax = mpc.bus(:,12);
voltMin = mpc.bus(:,13);

% return
result_pf_pinj = zeros(numOfBus,1);
result_pf_pinj([busNumOfSlack;busNumOfPV]) = result_pf.gen(:,2) - result_pf.bus([busNumOfSlack;busNumOfPV],3);
result_pf_pinj(busNumOfPQ) = - result_pf.bus(busNumOfPQ,3);

result_pf_qinj = zeros(numOfBus,1);
result_pf_qinj([busNumOfSlack;busNumOfPV]) = result_pf.gen(:,3) - result_pf.bus([busNumOfSlack;busNumOfPV],4);
result_pf_qinj(busNumOfPQ) = - result_pf.bus(busNumOfPQ,4);

result_pf_pinj = result_pf_pinj/mpc.baseMVA;
result_pf_qinj = result_pf_qinj/mpc.baseMVA;

%% MATLAB x = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon)

%% YALMIP 
% Parameters
c2 = mpc.gencost(:,5);
c1 = mpc.gencost(:,6);
c0 = mpc.gencost(:,7);
pd = mpc.bus(:,3);
qd = mpc.bus(:,4);


% Variables
e = sdpvar(numOfBus,1);
f = sdpvar(numOfBus,1);
pg = sdpvar(numOfBus,1);
qg = sdpvar(numOfBus,1);

% Constraints
constraint_yalmip = [];
constraint_yalmip = [constraint_yalmip, 0 <= qg(busNumOfSlack)  <= 10];
constraint_yalmip = [constraint_yalmip, f(busNumOfSlack) == 0];

for i = 1:numOfPVBuses
    busIndex = busNumOfPV(i);
    constraint_yalmip = [constraint_yalmip, 0 <= pg(busIndex) ];
end

for i = 1:numOfPQBuses
    busIndex = busNumOfPQ(i);
    constraint_yalmip = [constraint_yalmip, pg(busIndex) == 0];
    constraint_yalmip = [constraint_yalmip, qg(busIndex) == 0];
end

for i = 1:numOfBuses
    busIndex = i;
    constraint_yalmip = [constraint_yalmip, mpc.bus(busIndex,end)^2 <= e(busIndex)^2 + f(busIndex)^2 <= mpc.bus(busIndex,end-1)^2];
    constraint_yalmip = [constraint_yalmip, G(busIndex,:)*(e(busIndex)*e + f(busIndex)*f) + B(busIndex,:)*(f(busIndex)*e - e(busIndex)*f) == (pg(busIndex) - pd(busIndex))/mpc.baseMVA];
    constraint_yalmip = [constraint_yalmip, G(busIndex,:)*(f(busIndex)*e - e(busIndex)*f) - B(busIndex,:)*(e(busIndex)*e + f(busIndex)*f) == (qg(busIndex) - qd(busIndex))/mpc.baseMVA];
end

obj_yalmip = pg([busNumOfSlack;busNumOfPV])'*diag(c2)*pg([busNumOfSlack;busNumOfPV]) + c1'*pg([busNumOfSlack;busNumOfPV]) + sum(c0);

% % Assign power flow result as initial point
assign(e,result_pf.bus(:,8) .* cos(2*pi/360*result_pf.bus(:,9)));
assign(f,result_pf.bus(:,8) .* sin(2*pi/360*result_pf.bus(:,9)));
assign(pg([busNumOfSlack;busNumOfPV]),result_pf.gen(:,2));
assign(qg([busNumOfSlack;busNumOfPV]),result_pf.gen(:,3));

ops = sdpsettings('usex0',1);
optimize(constraint_yalmip,obj_yalmip,ops);

cost_yalmip = value(obj_yalmip);
[result_pf.gen(:,2), result_opf.gen(:,2), value(pg([busNumOfSlack;busNumOfPV])), result_pf.gen(:,3), result_opf.gen(:,3), value(qg([busNumOfSlack;busNumOfPV])), result_opf.gen(:,5), result_opf.gen(:,4)]
[mpc.bus(:,1), mpc.bus(:,2), value(e), value(f), sqrt(value(e).^2 + value(f).^2), mpc.bus(:,end), mpc.bus(:,end-1)]


%% Display Message
fprintf(' == MatPower Results ============================\n');
fprintf('       cost_total_pf = %f\n',cost_total_pf);
fprintf('      cost_total_opf = %f\n',cost_total_opf);
fprintf(' -- ACOPF (Minimum Cost) ------------------------\n');
fprintf('   cost_total_yalmip = %f\n',cost_yalmip);



fprintf(' ------------------------------------------------\n');
% fprintf(' ------------------------------------------------\n');
% fprintf(' ------------------------------------------------\n');
fprintf(' ================================================\n')