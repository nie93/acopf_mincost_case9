function [c,ceq] = acopfcon(x,mpc)
    numOfBuses = size(mpc.bus,1);
    numOfBranches = size(mpc.branch,1); 
    slackBusNum = find(mpc.bus(:,2) == 3);
    
    Ybus = makeYbus(mpc);
    G = real(Ybus);
    B = imag(Ybus);
    
    pg = x(1:numOfBuses);
    qg = x(numOfBuses+1:2*numOfBuses);
    e = x(2*numOfBuses+1:3*numOfBuses);
    f = x(3*numOfBuses+1:4*numOfBuses);
    
    pd = mpc.bus(:,3);
    qd = mpc.bus(:,4);    

    c = [];        
    ceq = [];
    
    % Imagine part of voltage on slack bus (f_{slack} == 0)
    ceq = [ceq; f(slackBusNum)];
    
    % Voltage Magnitudes
    for i = 1:numOfBuses
        c = [c; e(i)^2 + f(i)^2 - mpc.bus(i,12)^2];
        c = [c; mpc.bus(i,13)^2 - e(i)^2 + f(i)^2];
    end
    
    % Power flow balance
    for i = 1:numOfBuses
        ceq = [ceq; G(i,:)*(e(i)*e + f(i)*f) + B(i,:)*(f(i)*e - e(i)*f) - (pg(i) - pd(i))/mpc.baseMVA];
        ceq = [ceq; G(i,:)*(f(i)*e - e(i)*f) - B(i,:)*(e(i)*e + f(i)*f) - (qg(i) - qd(i))/mpc.baseMVA];
    end
    
    % Branch flow limits
    for i = 1:numOfBranches
        if (mpc.branch(i,6) ~= 0)
            fromBusIndex = mpc.branch(i,1);
            toBusIndex = mpc.branch(i,2);
            gij = G(fromBusIndex,toBusIndex);
            bij = B(fromBusIndex,toBusIndex);
            temp_mat = 0.5* [ -2*gij, gij, 0, -bij; gij, 0, bij, 0; 0, bij, -2*gij, gij; -bij, 0, gij, 0];
            c = [c; -mpc.branch(i,6) + mpc.baseMVA * abs([e(fromBusIndex); e(toBusIndex); f(fromBusIndex); f(toBusIndex)]' * temp_mat * [e(fromBusIndex); e(toBusIndex); f(fromBusIndex); f(toBusIndex)])];
            c = [c; -mpc.branch(i,6) + mpc.baseMVA * abs([e(toBusIndex); e(fromBusIndex); f(toBusIndex); f(fromBusIndex)]' * temp_mat * [e(toBusIndex); e(fromBusIndex); f(toBusIndex); f(fromBusIndex)])];
        end
    end
    
end