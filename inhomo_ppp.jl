function inhomogPPP(labmda,x1,x2,y1,y2,Nx,Ny;showplot=true)
    # Given an intensity function lambda, simulate an inhomogenous ppp

    lambda_max = maximum(lambda)
    p_homog = homogPPP(lambda_max,x1,x2,y1,y2,showplot=false)
    nmax = size(p_homog, 1)
    P = lambda/lambda_max;
    x = linspace(x1,x2,Nx)
    y = linspace(y1,y2,Ny)

    # Accpetance/rejection of points
    spl = Spline2D(x,y,P; kx=1, ky=1, s=0.0)
    p = zeros(nmax,2)
    cnt = 1
    for i=1:nmax
        p_keep = evaluate(spl, p_homog[i,1], p_homog[i,2])
        u = rand(1)[1]
        if u < p_keep
            p[cnt,:] = p_homog[i,:]
            cnt += 1
        end
    end
    p = p[1:(cnt-1),:]

    if showplot
        scatter(p[:,1],p[:,2],xlims=(x1,x2),ylims=(y1,y2))
    end

    return(p)
end

function inhomogPPPn(n,lambda,x1,x2,y1,y2,Nx,Ny;showplot=false)
    # Given an intensity function lambda, simulate an inhomogenous ppp with n points

    lambda_max = maximum(lambda)
    nmax = size(p_homog, 1)
    P = lambda/lambda_max;
    x = linspace(x1,x2,Nx)
    y = linspace(y1,y2,Ny)

    # Accpetance/rejection of points
    spl = Spline2D(x,y,P; kx=1, ky=1, s=0.0)
    p = zeros(n,2)
    np = 0
    while np < n
        x = uniform_pt_gen(1,x1,x2,y1,y2)
        p_keep = evaluate(spl, x[i,1], x[i,2])
        u = rand(1)[1]
        if u < p_keep
            p[cnt,:] = x
        end
    end

    if showplot
        scatter(p[:,1],p[:,2],xlims=(x1,x2),ylims=(y1,y2))
    end

    return(p)
end
