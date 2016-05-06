using Dierckx
using Distributions
using StatsFuns
using DataFrames
using Distances
using Plots
pyplot()
macro dlsym(func, lib)
    z, zlocal = gensym(string(func)), gensym()
    eval(current_module(),:(global $z = C_NULL))
    z = esc(z)
    quote
        let $zlocal::Ptr{Void} = $z::Ptr{Void}
            if $zlocal == C_NULL
               $zlocal = Libdl.dlsym($(esc(lib))::Ptr{Void}, $(esc(func)))
               global $z = $zlocal
            end
            $zlocal
        end
    end
end
