function cg = get_cost(mpc,pg)
    c2 = mpc.gencost(:,5);
    c1 = mpc.gencost(:,6);
    c0 = mpc.gencost(:,7);
    
    cg = pg'*diag(c2)*pg + c1'*pg + sum(c0);
end