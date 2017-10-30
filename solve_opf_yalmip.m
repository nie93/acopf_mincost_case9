function [yval, pg_opt, qg_opt]= solve_opf_yalmip(mpc)
    result_pf = runpf(mpc);
    result_opf = runopf(mpc);

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


    % Parameters
    c2 = mpc.gencost(:,5);
    c1 = mpc.gencost(:,6);
    c0 = mpc.gencost(:,7);
    pd = mpc.bus(:,3);
    qd = mpc.bus(:,4);


    % Variables
    e = sdpvar(numOfBuses,1);
    f = sdpvar(numOfBuses,1);
    pg = sdpvar(numOfBuses,1);
    qg = sdpvar(numOfBuses,1);

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

    for i = 1:numOfBranches
        fromBusIndex = mpc.branch(i,1);
        toBusIndex = mpc.branch(i,2);
        gij = G(fromBusIndex,toBusIndex);
        bij = B(fromBusIndex,toBusIndex);
        temp_mat = [ -2*gij, gij, 0, -bij; gij, 0, bij, 0; 0, bij, -2*gij, gij; -bij, 0, gij, 0];

        constraint_yalmip = [constraint_yalmip, mpc.branch(i,12) <= mpc.baseMVA * [e(fromBusIndex); e(toBusIndex); f(fromBusIndex); f(toBusIndex)]' * temp_mat * [e(fromBusIndex); e(toBusIndex); f(fromBusIndex); f(toBusIndex)] <= mpc.branch(i,13) ];
    end

    obj_yalmip = pg([busNumOfSlack;busNumOfPV])'*diag(c2)*pg([busNumOfSlack;busNumOfPV]) + c1'*pg([busNumOfSlack;busNumOfPV]) + sum(c0);

    % % Assign power flow result as initial point
    assign(e,result_pf.bus(:,8) .* cos(2*pi/360*result_pf.bus(:,9)));
    assign(f,result_pf.bus(:,8) .* sin(2*pi/360*result_pf.bus(:,9)));
    assign(pg([busNumOfSlack;busNumOfPV]),result_pf.gen(:,2));
    assign(qg([busNumOfSlack;busNumOfPV]),result_pf.gen(:,3));

    % % Assign flat start as initial point
    % assign(e,ones(numOfBus,1));
    % assign(f,zeros(numOfBus,1));
    % assign(pg([busNumOfSlack;busNumOfPV]),ones(numOfGen,1));
    % assign(qg([busNumOfSlack;busNumOfPV]),ones(numOfGen,1));

    ops = sdpsettings('usex0',1);
    optimize(constraint_yalmip,obj_yalmip,ops);

    cost_yalmip = value(obj_yalmip);
    [result_pf.gen(:,2), result_opf.gen(:,2), value(pg([busNumOfSlack;busNumOfPV])), ...
        result_pf.gen(:,3), result_opf.gen(:,3), value(qg([busNumOfSlack;busNumOfPV])), ...
        result_opf.gen(:,5), result_opf.gen(:,4)]

    [mpc.bus(:,1), mpc.bus(:,2), ...
        abs(value(e) + sqrt(-1)* value(f)), rad2deg(angle(value(e) + sqrt(-1)* value(f))), ...
        result_opf.bus(:,8), result_opf.bus(:,9), ...
        mpc.bus(:,end), mpc.bus(:,end-1)]

    lineflows = zeros(numOfBranches,2);

    for i = 1:numOfBranches
        fromBusIndex = mpc.branch(i,1);
        toBusIndex = mpc.branch(i,2);
        gij = G(fromBusIndex,toBusIndex);
        bij = B(fromBusIndex,toBusIndex);
        temp_mat = 0.5* [ -2*gij, gij, 0, -bij; gij, 0, bij, 0; 0, bij, -2*gij, gij; -bij, 0, gij, 0];
        ei = value(e(fromBusIndex));
        ej = value(e(toBusIndex));
        fi = value(f(fromBusIndex));
        fj = value(f(toBusIndex));
        lineflows(i,1) = [ei; ej; fi; fj]' * temp_mat * [ei; ej; fi; fj];
        lineflows(i,2) = [ej; ei; fj; fi]' * temp_mat * [ej; ei; fj; fi];
    end

    [mpc.branch(:,1), mpc.branch(:,2), mpc.branch(:,3), ...
        result_pf.branch(:,14), result_opf.branch(:,14), 100* lineflows, 100* abs(lineflows(:,1) + lineflows(:,2))]

    
    %% Outputs
    yval = cost_yalmip;
    pg_opt = value(pg([busNumOfSlack;busNumOfPV]));
    qg_opt = value(qg([busNumOfSlack;busNumOfPV]));
end