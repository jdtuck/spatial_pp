function thomasPP(mu, nu, sig; Ty=1., showplot=true)
    mu *= Ty
    d1 = Poisson(mu)
    d2 = Poisson(nu)
    d3 = Uniform()
    npts = rand(d1)

    xp = rand(d3,npts)
    yp = rand(d3,npts)*Ty
    ncl = zeros(Integer, npts)
    omax = quantile(d2,.9999)

    xcl = zeros(npts,omax)
    ycl = zeros(npts,omax)

    for i = 1:npts
        ncl[i] = round(Integer,rand(d2))
        for j = 1:ncl[i]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                r = sig*sqrt(-2*log(1-rand(d3)))
                theta = 2*pi*rand(d3)
                xclij = xp[i]+r*cos(theta)
                yclij = yp[i]+r*sin(theta)
            end

            xcl[i,j] = xclij
            ycl[i,j] = yclij
        end
    end

    xclf = Float64[]
    yclf = Float64[]
    for i=1:npts
        append!(xclf, vec(xcl[i,1:ncl[i]]))
        append!(yclf, vec(ycl[i,1:ncl[i]]))
    end

    if showplot
        scatter(xclf, yclf, xlims=(0,1), ylims=(0,Ty), ms=3)
        scatter!(xp, yp, mc=:red, ms=3)
    end

    return(xclf,yclf,xp,yp)
end

function inverse_powerPP(mu, nu, p, c; Ty=1., showplot=true)
    mu *= Ty
    d1 = Poisson(mu)
    d2 = Poisson(nu)
    d3 = Uniform()
    npts = rand(d1)

    xp = rand(d3,npts)
    yp = rand(d3,npts)*Ty
    ncl = zeros(Integer, npts)
    omax = quantile(d2,.9999)

    xcl = zeros(npts,omax)
    ycl = zeros(npts,omax)

    for ii = 1:npts
        ncl[ii] = rand(d2)
        for jj = 1:ncl[ii]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                r = c*((1-rand(d3))^(1/(1-p))-1)
                theta = 2*pi*rand(d3)
                xclij = xp[ii]+r*cos(theta)
                yclij = yp[ii]+r*sin(theta)
            end

            xcl[ii,jj] = xclij
            ycl[ii,jj] = yclij
        end
    end

    xclf = Float64[]
    yclf = Float64[]
    for ii=1:npts
        append!(xclf, vec(xcl[ii,1:ncl[ii]]))
        append!(yclf, vec(ycl[ii,1:ncl[ii]]))
    end

    if showplot
        scatter(xclf, yclf, xlims=(0,1), ylims=(0,Ty), ms=3)
        scatter!(xp, yp, mc=:red, ms=3)
    end

    return(xclf,yclf,xp,yp)
end

function thomasPP_genA(mu, nu, a, sigma_1, sigma_2; Ty=1,showplot=true)
    mu *= Ty
    d1 = Poisson(mu)
    d2 = Poisson(nu)
    d3 = Uniform()
    npts = rand(d1)

    xp = rand(d3,npts)
    yp = rand(d3,npts)*Ty
    ncl = zeros(Integer, npts)
    omax = quantile(d2,.9999)

    xcl = zeros(npts,omax)
    ycl = zeros(npts,omax)

    for i = 1:npts
        ncl[i] = rand(d2)
        for j = 1:ncl[i]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                U2 = rand(d3)
                U1 = rand(d3)
                if (U2<=a)
                    r = sigma_1*sqrt(-2*log(1-U1))
                else
                    r = sigma_2*sqrt(-2*log(1-U1))
                end
                theta = 2*pi*rand(d3)
                xclij = xp[i]+r*cos(theta)
                yclij = yp[i]+r*sin(theta)
            end

            xcl[i,j] = xclij
            ycl[i,j] = yclij
        end
    end

    xclf = Float64[]
    yclf = Float64[]
    for i=1:npts
        append!(xclf, vec(xcl[i,1:ncl[i]]))
        append!(yclf, vec(ycl[i,1:ncl[i]]))
    end

    if showplot
        scatter(xclf, yclf, xlims=(0,1), ylims=(0,Ty), ms=3)
        scatter!(xp, yp, mc=:red, ms=3)
    end

    return(xclf,yclf,xp,yp)
end

function thomasPP_genB(mu_1, mu_2, nu, sigma_1, sigma_2; Ty=1., showplot=true)
    d1 = Poisson(mu_1)
    d12 = Poisson(mu_2)
    d2 = Poisson(nu)
    d3 = Uniform()
    npts1 = rand(d1)
    npts2 = rand(d12)

    xp1 = rand(d3,npts1)
    yp1 = rand(d3,npts1)*Ty
    ncl1 = zeros(Integer, npts1)
    omax = quantile(d2,.9999)

    xcl1 = zeros(npts1,omax)
    ycl1 = zeros(npts1,omax)

    for i = 1:npts1
        ncl1[i] = rand(d2)
        for j = 1:ncl1[i]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                r = sigma_1*sqrt(-2*log(1-rand(d3)))
                theta = 2*pi*rand(d3)
                xclij = xp1[i]+r*cos(theta)
                yclij = yp1[i]+r*sin(theta)
            end

            xcl1[i,j] = xclij
            ycl1[i,j] = yclij
        end
    end

    xp2 = rand(d3,npts2)
    yp2 = rand(d3,npts2)*Ty
    ncl2 = zeros(Integer, npts2)
    omax = quantile(d2,.9999)

    xcl2 = zeros(npts2,omax)
    ycl2 = zeros(npts2,omax)

    for i = 1:npts2
        ncl2[i] = rand(d2)
        for j = 1:ncl2[i]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                r = sigma_2*sqrt(-2*log(1-rand(d3)))
                theta = 2*pi*rand(d3)
                xclij = xp2[i]+r*cos(theta)
                yclij = yp2[i]+r*sin(theta)
            end

            xcl2[i,j] = xclij
            ycl2[i,j] = yclij
        end
    end

    xclf = Float64[]
    yclf = Float64[]
    for i=1:npts1
        append!(xclf, vec(xcl1[i,1:ncl1[i]]))
        append!(yclf, vec(ycl1[i,1:ncl1[i]]))
    end
    for i=1:npts2
        append!(xclf, vec(xcl2[i,1:ncl2[i]]))
        append!(yclf, vec(ycl2[i,1:ncl2[i]]))
    end

    xp = [xp1;xp2]
    yp = [yp1;yp2]

    if showplot
        scatter(xclf, yclf, xlims=(0,1), ylims=(0,Ty), ms=3)
        scatter!(xp, yp, mc=:red, ms=3)
    end

    return(xclf,yclf,xp,yp)
end

function thomasPP_genC(mu_1, mu_2, nu_1, nu_2, sigma_1, sigma_2; Ty=1., showplot=true)
    d1 = Poisson(mu_1)
    d12 = Poisson(mu_2)
    d2 = Poisson(nu_1)
    d22 = Poisson(nu_2)
    d3 = Uniform()
    npts1 = rand(d1)
    npts2 = rand(d12)

    xp1 = rand(d3,npts1)
    yp1 = rand(d3,npts1)*Ty
    ncl1 = zeros(Integer, npts1)
    omax = quantile(d2,.9999)

    xcl1 = zeros(npts1,omax)
    ycl1 = zeros(npts1,omax)

    for i = 1:npts1
        ncl1[i] = rand(d2)
        for j = 1:ncl1[i]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                r = sigma_1*sqrt(-2*log(1-rand(d3)))
                theta = 2*pi*rand(d3)
                xclij = xp1[i]+r*cos(theta)
                yclij = yp1[i]+r*sin(theta)
            end

            xcl1[i,j] = xclij
            ycl1[i,j] = yclij
        end
    end

    xp2 = rand(d3,npts2)
    yp2 = rand(d3,npts2)*Ty
    ncl2 = zeros(Integer, npts2)
    omax = quantile(d22,.9999)

    xcl2 = zeros(npts2,omax)
    ycl2 = zeros(npts2,omax)

    for i = 1:npts2
        ncl2[i] = rand(d22)
        for j = 1:ncl2[i]
            xclij = -1
            yclij = -1
            while ((xclij<0 || xclij>1) || (yclij<0 || yclij>1))
                r = sigma_2*sqrt(-2*log(1-rand(d3)))
                theta = 2*pi*rand(d3)
                xclij = xp2[i]+r*cos(theta)
                yclij = yp2[i]+r*sin(theta)
            end

            xcl2[i,j] = xclij
            ycl2[i,j] = yclij
        end
    end

    xclf = Float64[]
    yclf = Float64[]
    for i=1:npts1
        append!(xclf, vec(xcl1[i,1:ncl1[i]]))
        append!(yclf, vec(ycl1[i,1:ncl1[i]]))
    end
    for i=1:npts2
        append!(xclf, vec(xcl2[i,1:ncl2[i]]))
        append!(yclf, vec(ycl2[i,1:ncl2[i]]))
    end

    xp = [xp1;xp2]
    yp = [yp1;yp2]

    if showplot
        scatter(xclf, yclf, xlims=(0,1), ylims=(0,Ty), ms=3)
        scatter!(xp, yp, mc=:red, ms=3)
    end

    return(xclf,yclf,xp,yp)
end

function thomas_palm(x, y, mu, nu, sigma, delta; Ty=1, showplot=true)
    np = length(x)
    R=0.5e0
    rmax=R/delta
    jmax = floor(Integer,rmax)
    nc = zeros(jmax)

    tx = 1.0
    ty = Ty
    rr, nn = bdry(x, y, tx, ty)

    for i=1:nn
        id = rr[i]/delta+1
        id = round(Integer,id)
        if (id <= jmax)
            nc[id] += 1
        end
    end
    delta1 = 0.001

    np_palm = zeros(jmax)
    palm_normal = zeros(jmax)
    for i=1:jmax
        r = delta1*i
        np_palm[i] = ((nc[i]/np)/((pi*(r+delta1)^2)-(pi*(r^2))))/np

        ae = exp(-r^2/(4*sigma^2))
        palm_normal[i] = (mu*nu+nu*ae/(4*pi*sigma^2))/(mu*nu)
    end

    r = collect(1:jmax)*delta
    if showplot
        ymax = max(maximum(np_palm),maximum(palm_normal))*2
        scatter(r,np_palm, xaxis=("r", :log),
                yaxis=("r", :log, (0.5,ymax)))
        plot!(r,palm_normal)
    end

    return(r, np_palm, palm_normal)
end

function inverse_power_palm(x, y, mu, nu, p, c, delta, x2; Ty=1., showplot=true)
    np = length(x)
    amu = zeros(1)
    anu = zeros(1)
    ap = zeros(1)
    ac = zeros(1)
    amu[1] = mu
    anu[1] = nu
    ap[1] = p
    ac[1] = c
    R = 0.5
    rmax = R/delta
    jmax = floor(Integer, rmax)
    np_palm = Array(Float64, jmax)
    palm_normal = Array(Float64, jmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:xqgausipf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
          Ptr{Float64}), x, y, &np, &delta, &Ty, &x2, amu, anu, ap, ac, &1,
          &jmax, np_palm, palm_normal)

    r = collect(1:jmax)*delta
    if showplot
        ymax = max(maximum(np_palm),maximum(palm_normal))*2
        scatter(r,np_palm, xaxis=("r", :log),
                yaxis=("r", :log, (0.5,ymax)))
        plot!(r,palm_normal)
    end

    return(r, np_palm, palm_normal)
end

function thomas_palm_genA(x, y, mu, nu, a, sigma_1, sigma_2, delta, x2; Ty=1.,
                          showplot=true)
    np = length(x)
    amu = zeros(1)
    anu = zeros(1)
    aa = zeros(1)
    asigma_1 = zeros(1)
    asigma_2 = zeros(1)
    amu[1] = mu
    anu[1] = nu
    aa[1] = a
    asigma_1[1] = sigma_1
    asigma_2[1] = sigma_2
    R = 0.5
    rmax = R/delta
    jmax = floor(Integer, rmax)
    np_palm = Array(Float64, jmax)
    palm_normal = Array(Float64, jmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:xqgausaf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Float64}), x, y, &np, &delta, &Ty, &x2, amu, anu, aa,
        asigma_1, asigma_2, &1, &jmax, np_palm, palm_normal)

    r = collect(1:jmax)*delta
    if showplot
      ymax = max(maximum(np_palm),maximum(palm_normal))*2
      scatter(r,np_palm, xaxis=("r", :log),
              yaxis=("r", :log, (0.5,ymax)))
      plot!(r,palm_normal)
    end

    return(r, np_palm, palm_normal)
end

function thomas_palm_genB(x, y, mu_1, mu_2, nu, sigma_1, sigma_2, delta; Ty=1.,
                          showplot=true)

    np = length(x)
    R=0.5e0
    rmax=R/delta
    jmax = floor(Integer,rmax)
    nc = zeros(jmax)

    tx = 1.0
    ty = Ty
    rr, nn = bdry(x, y, tx, ty)

    for i=1:nn
        id = rr[i]/delta+1
        id = round(Integer,id)
        if (id <= jmax)
            nc[id] += 1
        end
    end
    delta1 = 0.001

    lambda = nu*(mu_1+mu_2)
    a = mu_1/(mu_1+mu_2)

    np_palm = zeros(jmax)
    palm_normal = zeros(jmax)
    for i=1:jmax
        r = delta1*i
        np_palm[i] = ((nc[i]/np)/((pi*(r+delta1)^2)-(pi*(r^2))))/np

        ae1 = a/(sigma_1^2) * exp(-r^2/(4*sigma_1^2))
        ae2 = (1-a)/(sigma_2^2) * exp(-r^2/(4*sigma_2^2))

        palm_normal[i] = (lambda + (nu/(4*pi))*(ae1+ae2))/lambda
    end

    r = collect(1:jmax)*delta
    if showplot
        ymax = max(maximum(np_palm),maximum(palm_normal))*2
        scatter(r,np_palm, xaxis=("r", :log),
                yaxis=("r", :log, (0.5,ymax)))
        plot!(r,palm_normal)
    end

    return(r, np_palm, palm_normal)
end

function thomas_palm_genC(x, y, mu_1, mu_2, nu, sigma_1, sigma_2, delta; Ty=1.,
                          showplot=true)

    np = length(x)
    R=0.5e0
    rmax=R/delta
    jmax = floor(Integer,rmax)
    nc = zeros(jmax)

    tx = 1.0
    ty = Ty
    rr, nn = bdry(x, y, tx, ty)

    for i=1:nn
        id = rr[i]/delta+1
        id = round(Integer,id)
        if (id <= jmax)
            nc[id] += 1
        end
    end
    delta1 = 0.001

    lambda = mu_1*nu_1+mu_2*nu_2
    a = (mu_1*nu_1)/lambda

    np_palm = zeros(jmax)
    palm_normal = zeros(jmax)
    for i=1:jmax
        r = delta1*i
        np_palm[i] = ((nc[i]/np)/((pi*(r+delta1)^2)-(pi*(r^2))))/np

        ae1 = (a*nu_1)/(sigma_1^2) * exp(-r^2/(4*sigma_1^2))
        ae2 = ((1-a)*nu_2)/(sigma_2^2) * exp(-r^2/(4*sigma_2^2))

        palm_normal[i] = (lambda + (1/(4*pi))*(ae1+ae2))/lambda
    end

    r = collect(1:jmax)*delta
    if showplot
        ymax = max(maximum(np_palm),maximum(palm_normal))*2
        scatter(r,np_palm, xaxis=("r", :log),
                yaxis=("r", :log, (0.5,ymax)))
        plot!(r,palm_normal)
    end

    return(r, np_palm, palm_normal)
end

function thomas_mple(x, y, mu::Float64, nu::Float64, sigma::Float64; Ty=1., epsi=.1e-2, process=0,
                     showplot=true)
    # epsi - the optimization procedure is iterated at most 1000 times until stderr becomes smaller than epsi.
    # process
    # report the process of minimizing. Allowed values are
    # 0 :	no report.
    # 1 :	output the process of minimizing the negative Palm log-likelihood function until the values converge to the MPLE values for given data.
    # 2 :	output the process of optimizing by the simplex with the normalized parameters depending on pa. The actual estimates are obtained by the indicated x-values times pa.
    # 3 :	output the both processes

    np = length(x)
    ipflg = process
    if (process==0 && showplot)
        ipflg = 2
    end
    if (process==1 && showplot)
        ipflg = 3
    end
    n = 3
    itmax = 1000
    itmax1 = 1
    if (ipflg > 1)
        itmax1 = itmax+1
    end
    ipmax = itmax*2
    if (ipflg==0 || ipflg==2)
        ipmax = 1
    end

    fn = Array(Float64, ipmax)
    mple = Array(Float64, ipmax*n)
    xx = Array(Float64, itmax1*n)
    std = Array(Float64, itmax1)
    f = Array(Float64, itmax1)
    itr = Array(Int32, 1)
    nip = Array(Int32, 1)
    ipr = Array(Int32, ipmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:smplxthomf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}), x, y, &np, &Ty, &mu, &nu, &sigma, &epsi, &itmax,
        &itmax1, &ipmax, fn, mple, xx, std, f, itr, nip, ipr, &ipflg)

    nip = nip[1]
    ipri = ipr[1:nip]
    mple = reshape(mple,ipmax,n)
    if (nip != 1)
        mple = mple[1:nip,:]
    end
    mples = DataFrame(mu=mple[:,1], nu=mple[:,2], sigma=mple[:,3])

    it1 = itr[1]+1
    it2 = 1
    if (ipflg==2 || ipflg==3)
        it2 = it1
    end
    xx = reshape(xx,n,itmax1)

    if (process < 2)
        f = f[it2]
        para = xx[1:n,it2]
        std = std[it2]
    else
        f = f[1:it2]
        para = xx[1:n,1:it2]
        std = std[1:it2]
    end
    logLp = fn[1:nip]
    logLs = f
    stderr = std
    pa_normal = para
    pa_sim = [mu; nu; sigma].*pa_normal
    pa_sim = pa_sim'
    pa_sim = DataFrame(mu=pa_sim[:,1], nu=pa_sim[:,2], sigma=pa_sim[:,3])

    d_aic = AIC(-1*f[end],n)
    d_bic = BIC(-1*f[end],n,np)

    if showplot
      plot(logLp, title="Log-Likelihoood")
      plot!(logLs)
      plot(stderr, title="Standard Error")
    end

    out = pp_mple(logLp, mples, logLs, stderr, pa_normal, pa_sim, ipri, d_aic,
                  d_bic, "thomas")

    return(out)

end

function inverse_power_mple(x, y, mu::Float64, nu::Float64, p::Float64,
                            c::Float64, x2; Ty=1., skipi::Int64=1,
                            epsi::Float64=.1e-2, process=0, showplot=true)
    # epsi - the optimization procedure is iterated at most 1000 times until stderr becomes smaller than epsi.
    # process
    # report the process of minimizing. Allowed values are
    # 0 :	no report.
    # 1 :	output the process of minimizing the negative Palm log-likelihood function until the values converge to the MPLE values for given data.
    # 2 :	output the process of optimizing by the simplex with the normalized parameters depending on pa. The actual estimates are obtained by the indicated x-values times pa.
    # 3 :	output the both processes

    np = length(x)
    ipflg = process
    if (process==0 && showplot)
        ipflg = 2
    end
    if (process==1 && showplot)
        ipflg = 3
    end
    n = 4
    itmax = 1000
    itmax1 = 1
    if (ipflg > 1)
        itmax1 = itmax+1
    end
    ipmax = itmax*2
    if (ipflg==0 || ipflg==2)
        ipmax = 1
    end

    fn = Array(Float64, ipmax)
    mple = Array(Float64, ipmax*n)
    xx = Array(Float64, itmax1*n)
    std = Array(Float64, itmax1)
    f = Array(Float64, itmax1)
    itr = Array(Int32, 1)
    nip = Array(Int32, 1)
    ipr = Array(Int32, ipmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:smplxipf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}), x, y,
        &np, &skipi, &Ty, &mu, &nu, &p, &c, &x2, &epsi, &itmax, &itmax1, &ipmax,
        fn, mple, xx, std, f, itr, nip, ipr, &ipflg)

    nip = nip[1]
    ipri = ipr[1:nip]
    mple = reshape(mple,ipmax,n)
    if (nip != 1)
        mple = mple[1:nip,:]
    end
    mples = DataFrame(mu=mple[:,1], nu=mple[:,2], p=mple[:,3], c=mple[:,4])

    it1 = itr[1]+1
    it2 = 1
    if (ipflg==2 || ipflg==3)
        it2 = it1
    end
    xx = reshape(xx,n,itmax1)

    f = f[1:it2]
    para = xx[1:n,1:it2]
    std = std[1:it2]
    logLp = fn[1:nip]
    logLs = f
    stderr = std
    pa_normal = para
    pa_sim = [mu; nu; p; c].*pa_normal
    pa_sim = pa_sim'
    pa_sim = DataFrame(mu=pa_sim[:,1], nu=pa_sim[:,2], p=pa_sim[:,3], c=pa_sim[:,4])

    d_aic = AIC(-1*f[end],n)
    d_bic = BIC(-1*f[end],n,np)

    if showplot
        plot(logLp, title="Log-Likelihoood")
        plot!(logLs)
        plot(stderr, title="Standard Error")
    end

    out = pp_mple(logLp, mples, logLs, stderr, pa_normal, pa_sim, ipri, d_aic,
                  d_bic, "inverse_power")

    return(out)

end

function thomas_genA_mple(x, y, mu::Float64, nu::Float64, a::Float64,
                          sigma_1::Float64, sigma_2::Float64, x2; Ty=1.,
                          skipi::Int64=1, epsi::Float64=.1e-2, process=0,
                          showplot=true)
    # epsi - the optimization procedure is iterated at most 1000 times until stderr becomes smaller than epsi.
    # process
    # report the process of minimizing. Allowed values are
    # 0 :	no report.
    # 1 :	output the process of minimizing the negative Palm log-likelihood function until the values converge to the MPLE values for given data.
    # 2 :	output the process of optimizing by the simplex with the normalized parameters depending on pa. The actual estimates are obtained by the indicated x-values times pa.
    # 3 :	output the both processes

    np = length(x)
    ipflg = process
    if (process==0 && showplot)
        ipflg = 2
    end
    if (process==1 && showplot)
        ipflg = 3
    end
    n = 5
    itmax = 1000
    itmax1 = 1
    if (ipflg > 1)
        itmax1 = itmax+1
    end
    ipmax = itmax*2
    if (ipflg==0 || ipflg==2)
        ipmax = 1
    end

    fn = Array(Float64, ipmax)
    mple = Array(Float64, ipmax*n)
    xx = Array(Float64, itmax1*n)
    std = Array(Float64, itmax1)
    f = Array(Float64, itmax1)
    itr = Array(Int32, 1)
    nip = Array(Int32, 1)
    ipr = Array(Int32, ipmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:smplxaf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}), x, y, &np, &skipi, &Ty, &mu, &nu, &a, &sigma_1, &sigma_2,
        &x2, &epsi, &itmax, &itmax1, &ipmax, fn, mple, xx, std, f, itr, nip, ipr,
        &ipflg)

    nip = nip[1]
    ipri = ipr[1:nip]
    mple = reshape(mple,ipmax,n)
    if (nip != 1)
        mple = mple[1:nip,:]
    end
    mples = DataFrame(mu=mple[:,1], nu=mple[:,2], a=mple[:,3], sigma_1=mple[:,4],
                      sigma_2=mple[:,5])

    it1 = itr[1]+1
    it2 = 1
    if (ipflg==2 || ipflg==3)
        it2 = it1
    end
    xx = reshape(xx,n,itmax1)

    f = f[1:it2]
    para = xx[1:n,1:it2]
    std = std[1:it2]
    logLp = fn[1:nip]
    logLs = f
    stderr = std
    pa_normal = para
    pa_sim = [mu; nu; a; sigma_1; sigma_2].*pa_normal
    pa_sim = pa_sim'
    pa_sim = DataFrame(mu=pa_sim[:,1], nu=pa_sim[:,2], a=pa_sim[:,3],
                       sigma_1=pa_sim[:,4], sigma_2=pa_sim[:,5])

    d_aic = AIC(-1*f[end],n)
    d_bic = BIC(-1*f[end],n,np)

    if showplot
       plot(logLp, title="Log-Likelihoood")
       plot!(logLs)
       plot(stderr, title="Standard Error")
    end

    out = pp_mple(logLp, mples, logLs, stderr, pa_normal, pa_sim, ipri, d_aic,
                 d_bic, "thomas_genA")

    return(out)

end

function thomas_genB_mple(x, y, mu_1::Float64, mu_2::Float64, nu::Float64,
                          sigma_1::Float64, sigma_2::Float64; Ty=1.0,
                          epsi::Float64=.1e-2, process=0, showplot=true)
    # epsi - the optimization procedure is iterated at most 1000 times until stderr becomes smaller than epsi.
    # process
    # report the process of minimizing. Allowed values are
    # 0 :	no report.
    # 1 :	output the process of minimizing the negative Palm log-likelihood function until the values converge to the MPLE values for given data.
    # 2 :	output the process of optimizing by the simplex with the normalized parameters depending on pa. The actual estimates are obtained by the indicated x-values times pa.
    # 3 :	output the both processes

    np = length(x)
    ipflg = process
    if (process==0 && showplot)
        ipflg = 2
    end
    if (process==1 && showplot)
        ipflg = 3
    end
    n = 5
    itmax = 1000
    itmax1 = 1
    if (ipflg > 1)
        itmax1 = itmax+1
    end
    ipmax = itmax*2
    if (ipflg==0 || ipflg==2)
        ipmax = 1
    end

    fn = Array(Float64, ipmax)
    mple = Array(Float64, ipmax*n)
    xx = Array(Float64, itmax1*n)
    std = Array(Float64, itmax1)
    f = Array(Float64, itmax1)
    itr = Array(Int32, 1)
    nip = Array(Int32, 1)
    ipr = Array(Int32, ipmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:smplxbf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}), x, y, &np, &Ty, &mu_1,
        &mu_2, &nu, &sigma_1, &sigma_2, &epsi, &itmax, &itmax1, &ipmax, fn, mple,
        xx, std, f, itr, nip, ipr, &ipflg)

    nip = nip[1]
    ipri = ipr[1:nip]
    mple = reshape(mple,ipmax,n)
    if (nip != 1)
        mple = mple[1:nip,:]
    end
    mples = DataFrame(mu_1=mple[:,1], mu_2=mple[:,2], nu=mple[:,3],
                      sigma_1=mple[:,4], sigma_2=mple[:,5])

    it1 = itr[1]+1
    it2 = 1
    if (ipflg==2 || ipflg==3)
        it2 = it1
    end
    xx = reshape(xx,n,itmax1)

    f = f[1:it2]
    para = xx[1:n,1:it2]
    std = std[1:it2]
    logLp = fn[1:nip]
    logLs = f
    stderr = std
    pa_normal = para
    pa_sim = [mu_1; mu_2; nu; sigma_1; sigma_2].*pa_normal
    pa_sim = pa_sim'
    pa_sim = DataFrame(mu_1=pa_sim[:,1], mu_2=pa_sim[:,2], nu=pa_sim[:,3],
                       sigma_1=pa_sim[:,4], sigma_2=pa_sim[:,5])

    d_aic = AIC(-1*f[end],n)
    d_bic = BIC(-1*f[end],n,np)

    if showplot
      plot(logLp, title="Log-Likelihoood")
      plot!(logLs)
      plot(stderr, title="Standard Error")
    end

    out = pp_mple(logLp, mples, logLs, stderr, pa_normal, pa_sim, ipri, d_aic,
                  d_bic, "thomas_genB")

    return(out)

end

function thomas_genC_mple(x, y, mu_1::Float64, mu_2::Float64, nu_1::Float64,
                          nu_2::Float64, sigma_1::Float64, sigma_2::Float64;
                          Ty=1.0, epsi::Float64=.1e-2, process=0, showplot=true)
    # epsi - the optimization procedure is iterated at most 1000 times until stderr becomes smaller than epsi.
    # process
    # report the process of minimizing. Allowed values are
    # 0 :	no report.
    # 1 :	output the process of minimizing the negative Palm log-likelihood function until the values converge to the MPLE values for given data.
    # 2 :	output the process of optimizing by the simplex with the normalized parameters depending on pa. The actual estimates are obtained by the indicated x-values times pa.
    # 3 :	output the both processes

    np = length(x)
    ipflg = process
    if (process==0 && showplot)
        ipflg = 2
    end
    if (process==1 && showplot)
        ipflg = 3
    end
    n = 5
    itmax = 1000
    itmax1 = 1
    if (ipflg > 1)
        itmax1 = itmax+1
    end
    ipmax = itmax*2
    if (ipflg==0 || ipflg==2)
        ipmax = 1
    end

    fn = Array(Float64, ipmax)
    mple = Array(Float64, ipmax*n)
    xx = Array(Float64, itmax1*n)
    std = Array(Float64, itmax1)
    f = Array(Float64, itmax1)
    itr = Array(Int32, 1)
    nip = Array(Int32, 1)
    ipr = Array(Int32, ipmax)

    nscluster = Libdl.dlopen("deps/src/nscluster/libnscluster.dylib")
    ccall(@dlsym(:smplxcf_, nscluster), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}), x, y,
        &np, &Ty, &mu_1, &mu_2, &nu_1, &nu_2, &sigma_1, &sigma_2, &epsi, &itmax,
        &itmax1, &ipmax, fn, mple, xx, std, f, itr, nip, ipr, &ipflg)

    nip = nip[1]
    ipri = ipr[1:nip]
    mple = reshape(mple,ipmax,n)
    if (nip != 1)
        mple = mple[1:nip,:]
    end
    mples = DataFrame(lambda=mple[:,1], nu=mple[:,2], a=mple[:,3],
                      sigma_1=mple[:,4], sigma_2=mple[:,5])

    it1 = itr[1]+1
    it2 = 1
    if (ipflg==2 || ipflg==3)
        it2 = it1
    end
    xx = reshape(xx,n,itmax1)

    f = f[1:it2]
    para = xx[1:n,1:it2]
    std = std[1:it2]
    logLp = fn[1:nip]
    logLs = f
    stderr = std
    pa_normal = para
    lambda = mu_1*nu_1+mu_2*nu_2
    a = mu_1*nu_1/lambda
    pa_sim = [lambda; nu_1; a; sigma_1; sigma_2].*pa_normal
    pa_sim = pa_sim'
    pa_sim = DataFrame(lambda=pa_sim[:,1], nu=pa_sim[:,2], a=pa_sim[:,3],
                       sigma_1=pa_sim[:,4], sigma_2=pa_sim[:,5])

    d_aic = AIC(-1*f[end],n)
    d_bic = BIC(-1*f[end],n,np)

    if showplot
        plot(logLp, title="Log-Likelihoood")
        plot!(logLs)
        plot(stderr, title="Standard Error")
    end

    out = pp_mple(logLp, mples, logLs, stderr, pa_normal, pa_sim, ipri, d_aic,
                  d_bic, "thomas_genC")

    return(out)

end
