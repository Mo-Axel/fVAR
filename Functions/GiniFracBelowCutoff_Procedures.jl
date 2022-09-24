function ProbMassDiff_Cutoff(PhatDens, densvalues_ss, cutoff, xgrid)

    ngrid = length(xgrid);
    theta_sinh = 1;
    ygrid = 1/(2*theta_sinh)*(exp.(theta_sinh*xgrid) - exp.(-theta_sinh*xgrid));
    Jacobian = 1/2*(exp.(theta_sinh*xgrid) + exp.(-theta_sinh*xgrid));

    densvalues_diff = (PhatDens - densvalues_ss)./Jacobian;
    ygrid_diff      = ygrid[2:xn] - ygrid[1:xn-1];
    probmass_diff   = densvalues_diff[2:xn].*ygrid_diff;
    out = sum(probmass_diff[ygrid[2:xn].<cutoff]);

    return out

end


function GiniCoef(PhatDens, xgrid)

    ngrid = length(xgrid);
    theta_sinh = 1;
    ygrid = 1/(2*theta_sinh)*(exp.(theta_sinh*xgrid) - exp.(-theta_sinh*xgrid));
    Jacobian = 1/2*(exp.(theta_sinh*xgrid) + exp.(-theta_sinh*xgrid));

    densvalues_pmean   = PhatDens./Jacobian;
    # compute probability mass function
    ygrid_diff        = ygrid[2:ngrid] - ygrid[1:ngrid-1];
    probmassfcn       = densvalues_pmean[2:ngrid].*ygrid_diff;
    probmassfcn[1]    = probmassfcn[1]+densvalues_pmean[1]; # for mixed distribution

    # compute Gini coefficient
    Silag = 0;
    G     = 0;
    Si    = 1;
    for ii = 1:(ngrid-1)
        Si = sum(probmassfcn[1:ii].*ygrid[2:(1+ii)]);
        G  = G+probmassfcn[ii]*(Silag+Si);
        Silag = Si;
    end

    out = 1-G/Si;

return out

end
