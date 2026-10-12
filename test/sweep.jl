# test/sweep.jl — decode + build every SIF in the collection (no CUTEst
# comparison); writes a TSV of per-problem status for coverage bucketing.
using ExaModelsSIF
const DIR = expanduser("~/git/pkg/cutest-sif")
const OUT = expanduser("~/git/pkg/ExaModelsSIF.jl/sweep_results.tsv")

open(OUT, "w") do io
    for f in sort(readdir(DIR))
        endswith(f, ".SIF") || continue
        name = f[1:end-4]
        status, detail, nvar = "", "", -1
        P = nothing
        try
            P = sif_decode(joinpath(DIR, f))
            nvar = length(P.varnames)
        catch err
            status = "DECODE_FAIL"
            detail = sprint(showerror, err)
        end
        if isempty(status)
            if nvar > 200_000
                status = "SKIP_BUILD"
            else
                try
                    ExaModelsSIF.build_exa(P)
                    status = "OK"
                catch err
                    status = "BUILD_FAIL"
                    detail = sprint(showerror, err)
                end
            end
        end
        detail = replace(first(detail, 160), '\t' => ' ', '\n' => ' ')
        println(io, name, '\t', status, '\t', nvar, '\t', detail)
        flush(io)
    end
end
println("sweep done")
