# test/sweep.jl — decode + build every SIF in the collection (no CUTEst
# comparison); appends per-problem status to a TSV. Restartable: already-
# recorded problems are skipped, a crash leaves its name in current.txt and
# is recorded as CRASH on the next start. Exits with code 42 after CHUNK
# problems so the driver can restart Julia (bounds accumulated compilation
# memory on this shared host).
using ExaModelsSIF
const DIR = get(ENV, "SIF_DIR", "/home/sushin/.julia/artifacts/1fa01d4c3b344b8cd32b8ef78c0abeac9bb9cce2")  # canonical: CUTEst_jll artifact
const SHARD = parse(Int, get(ENV, "SWEEP_SHARD", "0"))
const NSHARDS = parse(Int, get(ENV, "SWEEP_NSHARDS", "1"))
const SUF = NSHARDS == 1 ? "" : ".$SHARD"
const OUT = expanduser("~/git/pkg/ExaModelsSIF.jl/sweep_results$(SUF).tsv")
const CUR = expanduser("~/git/pkg/ExaModelsSIF.jl/current$(SUF).txt")
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
        allf = sort(filter(f -> endswith(f, ".SIF"), readdir(DIR)))
    if haskey(ENV, "SWEEP_ONLY")
        want = Set(strip.(readlines(ENV["SWEEP_ONLY"])))
        allf = [f for f in allf if f[1:end-4] in want]
    end
        for (fi, f) in enumerate(allf)
            (fi - 1) % NSHARDS == SHARD || continue
            name = f[1:end-4]
            name in done && continue
            write(CUR, name)
            status, detail, nvar = "", "", -1
            P = nothing
            try
                P = sif_decode(joinpath(DIR, f))
                nvar = length(P.varnames)
            catch err
                if err isa ExaModelsSIF.ExternalFortranError
                    # hand-coded Fortran: out of scope for this package since
                    # 2026-10-14 (the oracle tier moved to ExaCUTEst.jl), so
                    # record the class rather than trying to evaluate it
                    status = "ORACLE_TIER"
                    detail = sprint(showerror, err)
                else
                    status = "DECODE_FAIL"
                    detail = sprint(showerror, err)
                end
            end
            if isempty(status)
                if nvar > 2_000_000
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
