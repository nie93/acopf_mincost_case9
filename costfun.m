function f = costfun(x,mpc)
    numOfBuses = size(mpc.bus,1);
    
    c2 = zeros(numOfBuses,1); 
    c1 = zeros(numOfBuses,1); 
    c0 = zeros(numOfBuses,1); 
    
    c2(mpc.bus(:,2) ~= 1,1) = mpc.gencost(:,5);
    c1(mpc.bus(:,2) ~= 1,1) = mpc.gencost(:,6);
    c0(mpc.bus(:,2) ~= 1,1) = mpc.gencost(:,7);
    
    pg = x(1:numOfBuses);
    f = pg'*diag(c2)*pg + c1'*pg + sum(c0);
    
end