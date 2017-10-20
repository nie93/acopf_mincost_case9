function cg = get_cost_pinj(mpc,pinj)
    numOfBus= size(mpc.bus,1);
    pinj = pinj*mpc.baseMVA;
    pd = mpc.bus(:,3);
    c2 = zeros(numOfBus,1);
    c1 = zeros(numOfBus,1);
    c0 = zeros(numOfBus,1);
    
    busNumOfSlack = find(mpc.bus(:,2) == 3);
    busNumOfPV = find(mpc.bus(:,2) == 2);
    busNumOfPQ = find(mpc.bus(:,2) == 1);
    
    c2([busNumOfSlack;busNumOfPV]) = mpc.gencost(:,5);
    c1([busNumOfSlack;busNumOfPV]) = mpc.gencost(:,6);
    c0([busNumOfSlack;busNumOfPV]) = mpc.gencost(:,7);
    
    c2(busNumOfPQ) = 0;
    c1(busNumOfPQ) = 0;
    c0(busNumOfPQ) = 0;
    
    pg = (pinj+pd);
    
    cg = pg'*diag(c2)*pg + c1'*pg + sum(c0);
end