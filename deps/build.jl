@unix_only begin
    cd(joinpath(dirname(@__FILE__), "src", "nscluster"))
    suffix = @osx? "dylib" : "so"
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
end

@windows_only begin
    # these binaries were cross-compiled from Cygwin for x86_64 only using
    # the Makefile_win in the corresponding src directories and the windows bin
    # TODO: Build
    #run(`curl -LO https://github.com/jdtuck/spatial_ppp.jl/releases/download/v0.1.0/nscluster.7z`)
    #run(`7z x -y nscluster.7z`)
end
