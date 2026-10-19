# test/valsweep.jl — full-collection VALIDATION sweep: every problem's values
# and derivatives against the Fortran reference via test/pilot.jl's check().
# Sharded/restartable like test/sweep.jl. External-Fortran problems cannot be
# decoded and are recorded as EXTERNAL_FORTRAN without comparison.
#
#   julia --project=test test/valsweep.jl
include(joinpath(@__DIR__, "pilot.jl"))

const DIR2 = ENV["MASTSIF"]
const SHARD = parse(Int, get(ENV, "SWEEP_SHARD", "0"))
const NSHARDS = parse(Int, get(ENV, "SWEEP_NSHARDS", "1"))
const SUF = ".val$SHARD"
const OUT = joinpath(@__DIR__, "..", "val_results$(SUF).tsv")
const CUR = joinpath(@__DIR__, "..", "current$(SUF).txt")
const CHUNK = 60

done = Set{String}()
isfile(OUT) && for ln in eachline(OUT)
    push!(done, first(split(ln, '\t')))
end
if isfile(CUR)
    prev = strip(read(CUR, String))
    if !isempty(prev) && !(prev in done)
        open(OUT, "a") do io
            println(io, prev, "\tCRASH\tprocess died during validation")
        end
        push!(done, String(prev))
    end
end

function run_valsweep()
    nproc = 0
    open(OUT, "a") do io
        allf = sort(filter(f -> endswith(f, ".SIF"), readdir(DIR2)))
        if haskey(ENV, "SWEEP_ONLY")
            want = Set(strip.(readlines(ENV["SWEEP_ONLY"])))
            allf = [f for f in allf if f[1:end-4] in want]
        end
        for (fi, f) in enumerate(allf)
            (fi - 1) % NSHARDS == SHARD || continue
            name = f[1:end-4]
            name in done && continue
            write(CUR, name)
            status, detail = "", ""
            try
                sif_decode(joinpath(DIR2, f))   # cheap external-tier probe
                msgs = check(name)
                if isempty(msgs)
                    status = "PASS"
                else
                    status = "FAIL"
                    detail = join(msgs[1:min(end, 3)], " | ")
                end
            catch err
                if err isa ExaModelsSIF.ExternalFortranError
                    status = "EXTERNAL_FORTRAN"
                else
                    status = "ERROR"
                    detail = sprint(showerror, err)
                end
            end
            detail = replace(first(detail, 200), '\t' => ' ', '\n' => ' ')
            println(io, name, '\t', status, '\t', detail)
            flush(io)
            write(CUR, "")
            nproc += 1
            if nproc >= CHUNK
                println("chunk done")
                exit(42)
            end
        end
    end
    println("valsweep done")
end

run_valsweep()
