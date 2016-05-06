type pp_mple
    logLp::Array{Float64,1}
    mples::DataFrames.DataFrame
    logLs::Array{Float64,1}
    stderr::Array{Float64,1}
    pa_normal::Array{Float64,2}
    pa_sim::DataFrames.DataFrame
    ipri::Array{Int32,1}
    AIC::Float64
    BIC::Float64
    model::ASCIIString
end

function circle(sx,sy,cx,cy,r)
    x,y = meshgrid(-(cx-1):(sx-cx), -(cy-1):(sy-cy))
    c_mask = ((x.^2+y.^2)<=r^2)
    return(c_mask)
end

function bdry(x,y,tx,ty)
    np = length(x)
    rr = zeros(np*np)
    t1 = 0.5
    nn = 0
    for i = 1:np
        for j = 1:np
            if (j==i)
                continue
            end
            xx = x[j] - x[i]
            if (xx > tx/2)
                xx = -(tx-xx)
            end
            if (xx < -tx/2)
                xx = tx+xx
            end
            yy = y[j] - y[i]
            if (yy > ty/2)
                yy = -(ty-yy)
            end
            if (yy < -ty/2)
                yy = ty+yy
            end
            if (abs(xx)>t1 || abs(yy)>t1)
                continue
            end
            r2 = sqrt(xx^2+yy^2)
            if (r2 > t1)
                continue
            end
            nn += 1
            rr[nn] = r2
        end
    end

    return(rr,nn)
end

function calc_NN_dists(p)
    R = pairwise(Euclidean(), p')
    R_sort = sort(R,2)
    d = R_sort(:,2)
    return(d)
end

function birth_death_move(C,lambda,Nx,Ny,x1,x2,y1,y2)
    # Perform a birth, a death, or a move to create a candidate

    N = size(C,1)
    u = rand(1)[1]
    C_cand = C
    if N == 1
        if u<(1/2) # birth
            p = inhomogPPPn(1,lambda,x1,x2,y1,y2,Nx,Ny)
            C_cand = hcat(C_cand,p)
        else  # move
            idx = rand(1:N)
            ci = C_cand[idx,:]
            ci_xind = round((ci[1]-x1)*(Nx-1))+1
            ci_yind = round((ci[2]-y1)*(Ny-1))+1
            r = 0.1*minimum([Nx,Ny])
            B_ci = circle(Nx,Ny,ci_yind,ci_xind,r)
            lambda_ci = lambda.*B_ci'
            p = inhomogPPPn(1,lambda_ci,x1,x2,y1,y2,Nx,Ny)
            C_cand[idx,:] = p
        end
    else
        if u<(1/3)  # birth
            p = inhomogPPPn(1,lambda,x1,x2,y1,y2,Nx,Ny)
            C_cand = hcat(C_cand,p)
        elseif (u>(1/3) && u<(2/3))  # death
            idx = rand(1:N)
            C_cand = hcat(C_cand[1:(idx-1),:], C_cand[(idx+1):end,:])
        else  # move
            idx = rand(1:N)
            ci = C_cand[idx,:]
            ci_xind = round((ci[1]-x1)*(Nx-1))+1
            ci_yind = round((ci[2]-y1)*(Ny-1))+1
            r = 0.1*minimum([Nx,Ny])
            B_ci = circle(Nx,Ny,ci_yind,ci_xind,r)
            lambda_ci = lambda.*B_ci'
            p = inhomogPPPn(1,lambda_ci,x1,x2,y1,y2,Nx,Ny)
            C_cand[idx,:] = p
        end
    end

    return(C_cand)
end

function lambda_est(p,c,x1,x2,y1,y2,Nx,Ny,h;showplot=true)
    # This function computes a nonparametric, kernel-based estimation of the
    # intensity function of a point process p in a rectangular window. The
    # estimation uses a quartic kernel function defined over a circular domain.
    #
    # Inputs
    #   p:  An n x 2 matrix consisting of the coordinates of n points in 2D,
    #   c:  An n x 1 vector of edge correction coefficients. This vector is
    #       computed using the function edgeCorr.m. If no edge correction is
    #       needed or desired, input c=ones(n,1).
    #   x1, x2, y1, y2: The x and y coordinates of the boundaries of the
    #       rectangular viewing window [x1,x2] x [y1,y2]. These values may be
    #       in feet, meters, miles, etc.
    #   Nx, Ny: The number of x and y grid points where the intensity function
    #       is computed
    #   h:  The kernel bandwidth, i.e. the width of the circular domain on which
    #       it is defined. This value controls how smooth or rough the
    #       intensity function estimate is. A larger value of h leads to a
    #       smoother intensity function but at a cost of losing small features.
    #       For a square viewing window, consider values of h that are 5 to 15
    #       percent the side length.
    #   showplot: show plots
    #
    # Output
    #   lambda_hat: An Nx x Ny matrix that represents the values at each
    #       gridpoint of the estimated intensity function in the viewing
    #       window [x1,x2] x [y1,y2].
    n = size(p,1)
    x = linspace(x1,x2,Nx)
    y = linspace(y1,y2,Ny)
    lambda_hat = zeros(Nx,Ny)


    for i=1:Nx
        for j=1:Ny
            gridpt = [x[i], y[j]]
            d = zeros(n)
            for k=1:n
                d[k] = norm(p[k,:]-gridpt)
            end

            idx = find(d.<h)
            d2 = d[d.<h].^2
            if !isempty(d2)
                for ii = 1:length(idx)
                    m = idx[ii]
                    lambda_hat[i,j] = lambda_hat[i,j]+((3/(pi*h^2))*(1-d2[ii]/(h^2))^2)/c[m]
                end
            end
        end
    end

    if showplot
        contour(x,y,lambda_hat',fill=true,xlims=(x1,x2),ylims=(y1,y2))
    end

    return(labmda_hat)
end

function CSR_test(p; Nx=10, Ny=10)
    # this function tests for complete spatial randomness (CSR) of a given spatial point process
    # Nx and Ny partion W
    k = Nx*Ny
    qcount = zeros(Nx,Ny)

    # place each point in a quadrant and add to qcount
    for i=1:n
        ptrans = (p[i,:]+1)/2.
        xind = ceil(Integer, ptrans[1]*Nx)
        yind = ceil(Integer, ptrans[2]*Ny)
        qcount[xind,yind] += 1
    end

    # Calculate and test the index-of-dispersion
    ppq = reshape(qcount,k,1)
    I = (k-1)*var(ppq)/mean(ppq)
    c1 = chisqinvcdf(k-1,0.05)
    c2 = chisqinvcdf(k-1,0.95)
    reject = (I<c1)||(I>c2)

    return(I, reject)
end

function AIC(L, k)
    d = -2*log(L)+2*k
    return(d)
end

function BIC(L, k, n)
    d = -2*log(L)+k*log(n)
    return(d)
end
