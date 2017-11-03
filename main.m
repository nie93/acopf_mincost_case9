clear; clc;

%% Get data from MatPower
% mpc = loadcase('case9');
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
numOfBranches = size(mpc.branch,1); 
voltMax = mpc.bus(:,12);
voltMin = mpc.bus(:,13);

% lineRatings = [100 50 100 35 50 50 50 50 50 50 70 50 50 50 35 30 50 70 30 30]';
% mpc = set_linelimits(mpc,lineRatings);

% return
result_pf_pinj = zeros(numOfBuses,1);
result_pf_pinj([busNumOfSlack;busNumOfPV]) = result_pf.gen(:,2) - result_pf.bus([busNumOfSlack;busNumOfPV],3);
result_pf_pinj(busNumOfPQ) = - result_pf.bus(busNumOfPQ,3);

result_pf_qinj = zeros(numOfBuses,1);
result_pf_qinj([busNumOfSlack;busNumOfPV]) = result_pf.gen(:,3) - result_pf.bus([busNumOfSlack;busNumOfPV],4);
result_pf_qinj(busNumOfPQ) = - result_pf.bus(busNumOfPQ,4);

result_pf_pinj = result_pf_pinj/mpc.baseMVA;
result_pf_qinj = result_pf_qinj/mpc.baseMVA;

%% MATLAB x = fmincon(fun,x0,[],[],[],[],lb,ub,nonlcon)
lb_p = zeros(numOfBuses,1);
ub_p = zeros(numOfBuses,1);
lb_q = zeros(numOfBuses,1);
ub_q = zeros(numOfBuses,1);

lb_p(mpc.bus(:,2)~=1,1) = mpc.gen(:,10);
ub_p(mpc.bus(:,2)~=1,1) = mpc.gen(:,9);
lb_q(mpc.bus(:,2)~=1,1) = mpc.gen(:,5);
ub_q(mpc.bus(:,2)~=1,1) = mpc.gen(:,4);

lb = [lb_p; lb_q];
ub = [ub_p; ub_q];

x0 = zeros(4*numOfBuses,1);
x0(mpc.bus(:,2)~=1,1) = mpc.gen(:,2);
x0(find(mpc.bus(:,2)~=1) + numOfBuses,1) = mpc.gen(:,3);
x0(2*numOfBuses+1:3*numOfBuses) = 1;

options = optimoptions('fmincon','Display','off');
[x_opf,cost_opf,exitflag,output] = fmincon(@(x)costfun(x,mpc),x0,[],[],[],[],lb,ub,@(x)acopfcon(x,mpc),options);

pg_opf = x_opf(1:numOfBuses);
qg_opf = x_opf(numOfBuses+1:2*numOfBuses);
e_opf = x_opf(2*numOfBuses+1:3*numOfBuses);
f_opf = x_opf(3*numOfBuses+1:4*numOfBuses);

v_opf = abs(e_opf + sqrt(-1) * f_opf);
delta_opf = angle(e_opf + sqrt(-1) * f_opf);
deltadeg_opf = rad2deg(delta_opf);


%% Display
output
x_opf(mpc.bus(:,2)~=1,1)
[(1:14)', pg_opf, lb_p, ub_p, qg_opf, lb_q, ub_q, e_opf, f_opf, v_opf, deltadeg_opf]
[[(1:numOfBuses)';(1:numOfBuses)';(1:numOfBuses)';(1:numOfBuses)'],x_opf,x0]

lineflows = zeros(numOfBranches,2);

for i = 1:numOfBranches
    fromBusIndex = mpc.branch(i,1);
    toBusIndex = mpc.branch(i,2);
    gij = G(fromBusIndex,toBusIndex);
    bij = B(fromBusIndex,toBusIndex);
    temp_mat = 0.5* [ -2*gij, gij, 0, -bij; gij, 0, bij, 0; 0, bij, -2*gij, gij; -bij, 0, gij, 0];
    ei = e_opf(fromBusIndex);
    ej = e_opf(toBusIndex);
    fi = f_opf(fromBusIndex);
    fj = f_opf(toBusIndex);
    lineflows(i,1) = [ei; ej; fi; fj]' * temp_mat * [ei; ej; fi; fj];
    lineflows(i,2) = [ej; ei; fj; fi]' * temp_mat * [ej; ei; fj; fi];
end

% [mpc.branch(:,1), mpc.branch(:,2), mpc.branch(:,3), ...
%     result_pf.branch(:,14), result_opf.branch(:,14), ...
%     mpc.baseMVA * lineflows(:,1), mpc.baseMVA * lineflows(:,2), ...
%     lineRatings]

%% YALMIP 
% [cost_yalmip, pg, qg] = solve_opf_yalmip(mpc);

%% Display Message
fprintf(' == MatPower Results ============================\n');
fprintf('       cost_total_pf = %f\n',cost_total_pf);
fprintf('      cost_total_opf = %f\n',cost_total_opf);
fprintf(' -- ACOPF (Minimum Cost) fmincon ----------------\n');
fprintf('  cost_total_fmincon = %f\n',cost_opf);
fprintf(' -- ACOPF (Minimum Cost) YALMIP -----------------\n');
% fprintf('   cost_total_yalmip = %f\n',cost_yalmip);


fprintf(' ------------------------------------------------\n');
% fprintf(' ------------------------------------------------\n');
% fprintf(' ------------------------------------------------\n');
fprintf(' ================================================\n')