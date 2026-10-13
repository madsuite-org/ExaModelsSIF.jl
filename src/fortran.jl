# fortran.jl — translate the restricted Fortran expression language used in the
# ELEMENTS/GROUPS INDIVIDUALS sections of a SIF file into Julia ASTs, and
# evaluate those ASTs at ExaModels trace time over symbolic nodes.

"Parse a SIF/Fortran expression string into a Julia AST."
function fparse(s::AbstractString)
    t = replace(String(s), "**" => "^")
    # Fortran D exponents: 1.0D0, 2.D-3, .5D+2 -> e
    t = replace(t, r"(?<=[0-9.])[dD](?=[+-]?[0-9])" => "e")
    ex = Meta.parse(t)
    return ex
end

# Fortran intrinsic name (uppercased) -> Julia function
const JLFUN = Dict{Symbol,Function}(
    :EXP => exp, :LOG => log, :LOG10 => log10,
    :SIN => sin, :COS => cos, :TAN => tan,
    :ASIN => asin, :ACOS => acos, :ATAN => atan, :ATAN2 => atan,
    :SINH => sinh, :COSH => cosh, :TANH => tanh,
    :SQRT => sqrt, :ABS => abs,
    :MAX => max, :MIN => min, :MOD => mod,
    :ARCTAN => atan, :ARCSIN => asin, :ARCCOS => acos,
)

"""
    feval(ast, env)

Evaluate a parsed Fortran expression over `env::Dict{Symbol,Any}` mapping
symbol names to numbers or ExaModels nodes. Arithmetic operators and the
intrinsics in `JLFUN` are applied directly, so the result is an ExaModels
expression tree when the env contains nodes.
"""
function feval(e, env::AbstractDict{Symbol,<:Any})
    if e isa Symbol
        haskey(env, e) && return env[e]
        u = Symbol(uppercase(String(e)))
        for (k, v) in env
            Symbol(uppercase(String(k))) === u && return v
        end
        error("SIF expression references undefined symbol $e")
    elseif e isa Number
        return e isa Integer ? e : Float64(e)
    elseif e isa Expr
        if e.head === :call
            op = e.args[1]
            args = Any[feval(a, env) for a in e.args[2:end]]
            return fapply(op, args)
        else
            error("unsupported expression head $(e.head) in SIF expression")
        end
    else
        error("unsupported token $e in SIF expression")
    end
end

function fapply(op::Symbol, args)
    if op === :+
        return length(args) == 1 ? args[1] : reduce(+, args)
    elseif op === :-
        return length(args) == 1 ? -args[1] : args[1] - args[2]
    elseif op === :*
        return reduce(*, args)
    elseif op === :/
        return args[1] / args[2]
    elseif op === :^
        p = args[2]
        if p isa Number && isinteger(p)
            return args[1]^Int(p)
        end
        return args[1]^p
    end
    f = get(JLFUN, Symbol(uppercase(String(op))), nothing)
    f === nothing && error("unsupported function $op in SIF expression")
    return f(args...)
end
