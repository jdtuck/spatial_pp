function uniform_pt_gen(n, x1, x2, y1, y2)
    # Generates n points uniformly on W = [x1,x2] x [y1, y2]
    p = zeros(n,2)
    p[:,1] = (x2-x1)*rand(n)+x1
    p[:,2] = (y2-y1)*rand(n)+y1
    return(p)
end

function homogPPP(lambda,x1,x2,y1,y2;showplot=true)
    # Simulate a homogenous Poisson process on W = [x1,x2] x [y1,y2]
    area = (x2-x1)*(y2-y1)
    d = Poisson(lambda*area)
    n = rand(d)
    p = uniform_pt_gen(n,x1,x2,y1,y2)

    if showplot
        scatter(p[:,1],p[:,2], xlims=(x1,x2), ylims=(y1,y2))
    end
    return(p)
end

function homogPPP_est(p,x1,x2,y1,y2)
    n = size(p,1)
    area = (x2-x1)*(y2-y1)
    lambda_hat = n/area
    CI=[(1.96/2-sqrt(n))^2/area; (1.96/2+sqrt(n+1))^2/area]
    return(lambda_hat, CI)
end

function homogPPP_estNN(p)
    # calculate all nearest neighbor distances, d
    d = calc_NN_dists(p)

    # estimate lambda
    r_hat = mean(d.^2)
    lambda_hat = 1/(pi*r_hat)

    return(lambda_hat)
end
