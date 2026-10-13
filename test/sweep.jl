# test/sweep.jl — decode + build every SIF in the collection (no CUTEst
# comparison); appends per-problem status to a TSV. Restartable: already-
# recorded problems are skipped, a crash leaves its name in current.txt and
# is recorded as CRASH on the next start. Exits with code 42 after CHUNK
# problems so the driver can restart Julia (bounds accumulated compilation
# memory on this shared host).
using ExaModelsSIF
const DIR = expanduser("~/git/pkg/cutest-sif")
const OUT = expanduser("~/git/pkg/ExaModelsSIF.jl/sweep_results.tsv")
const CUR = expanduser("~/git/pkg/ExaModelsSIF.jl/current.txt")
const CHUNK = 150

done = Set{String}()
isfile(OUT) && for ln in eachline(OUT)
    push!(done, first(split(ln, '\t')))
end
if isfile(CUR)
    prev = strip(read(CUR, String))
    if !isempty(prev) && !(prev in done)
        open(OUT, "a") do io
            println(io, prev, "\tCRASH\t-1\tprocess died while processing this problem")
        end
        push!(done, String(prev))
    end
end

function run_sweep()
    nproc = 0
    open(OUT, "a") do io
        for f in sort(filter(f -> endswith(f, ".SIF"), readdir(DIR)))
            name = f[1:end-4]
            name in done && continue
            write(CUR, name)
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
            write(CUR, "")
            nproc += 1
            if nproc >= CHUNK
                println("chunk done")
                exit(42)
            end
        end
    end
    println("sweep done")

end

run_sweep()
