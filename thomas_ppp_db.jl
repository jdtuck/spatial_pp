function thomasPP(alpha,kappa,omega,x1,x2,y1,y2;showplot=true)
    # Simulate parents
    C = homogPPP(kappa,x1,x2,y1,y2,showplot=false)
    nC = size(C,1)

    # Create intesity function
    Nx = 101
    Ny = 101
    x = linspace(x1,x2,Nx)
    y = linspace(y1,y2,Ny)
    X,Y = meshgrid(x,y)
    xx = reshape(X,Nx*Ny,1)
    yy = reshape(Y,Nx*Ny,1)
    gridpoints = hcat(xx,yy)
    lambda = zeros(Nx,Ny)
    sigma = omega^2*eye(2)
    for i=1:nC
        x0 = C[i,:]
        d = MvNormal(x0,sigma)
        rho = pdf(d,gridpoints')
        rho = reshape(rho',Nx,Ny)
        rho *= alpha
        labmda += rho
    end

    # simulate children
    X = inhomogPPP(lambda,x1,x2,y1,y2,Nx,Ny,showplot=false)

    if showplot
        contour(x,y,lambda,fill=true)
        scatter(X[:,1],X[:,2],xlims=(x1,x2),ylims=(y1,y2), grid=true)
        scatter!(C[:,1], C[:,2], markercolor=:red)
    end

    return(X,C,lambda)
end

function thomasPP_inhomog(alpha,kappa,omega,beta,x1,x2,y1,y2;showplot=true, dx=0.01)
    # simulates a homogenous thomas process with inhomogenoues background clutter

    # simulate homogenous thomas process
    X1,C,lambda = thomasPP(alpha,kappa,omeaga,x1,x2,y1,y2,showplot=false)

    # define spatial covariates for background clutter
    x = x1:dx:x2
    y = y1:dx:y2
    Nx = length(x)
    Ny = length(y)
    XX, YY = meshgrid(x,y)
    s = zeros(Nx,Ny)
    s += beta[1]*ones(Nx,Ny)
    s += beta[2]*XX
    s += beta[3]*YY
    eta = exp(s)

    # sample from eta to obtain background clutter
    X2 = inhomogPPP(eta,x1,x2,y1,y2,Nx,Ny,showplot=false)

    # form entire child process
    X = vcat(X1,X2)
    nX = size(X,2)

    if showplot
        contour(x,y,eta,fill=true)
        scatter(X1[:,1],X1[:,2],xlims=(x1,x2),ylims=(y1,y2), grid=true)
        scatter!(X2[:,1],X2[:,2], makercolor=:red)
        scatter!(C[:,1],C[:,2], markercolor=:green)
    end

    return(X,C,lambda)
end
