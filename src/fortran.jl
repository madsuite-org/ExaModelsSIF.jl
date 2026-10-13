# fortran.jl — translate the restricted Fortran expression language used in the
# ELEMENTS/GROUPS INDIVIDUALS sections of a SIF file into Julia ASTs, and
# evaluate those ASTs at ExaModels trace time over symbolic nodes.

"Parse a SIF/Fortran expression string into a Julia AST."
function fparse(s::AbstractString)
    t = replace(String(s), "**" => "^")
    # Fortran logical operators (used in conditional temporaries)
    t = replace(t,
        r"\.LE\."i => " <= ", r"\.LT\."i => " < ",
        r"\.GE\."i => " >= ", r"\.GT\."i => " > ",
        r"\.EQ\."i => " == ", r"\.NE\."i => " != ",
        r"\.AND\."i => " && ", r"\.OR\."i => " || ",
        r"\.NOT\."i => " !", r"\.TRUE\."i => " true ", r"\.FALSE\."i => " false ")
    # Fortran D exponents: 1.0D0, 2.D-3, .5D+2 -> e
    t = replace(t, r"(?<=[0-9.])[dD](?=[+-]?[0-9])" => "e")
    # Fortran allows blanks between a function name and its argument list
    t = replace(t, r"([A-Za-z][A-Za-z0-9_]*)\s+\(" => s"\1(")
    # ... and inside numeric literals ("1. 0" is 1.0)
    t = replace(t, r"([0-9]\.)\s+([0-9])" => s"\1\2")
    t = replace(t, r"([0-9])\s+(\.[0-9])" => s"\1\2")
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
    :DBLE => identity, :FLOAT => identity, :REAL => identity,
    :SIGN => (a, b) -> abs(a) * (b < 0 ? -1.0 : 1.0),
    :INT => x -> trunc(x), :AINT => x -> trunc(x), :NINT => x -> round(x),
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
    elseif e isa Bool
        return e
    elseif e isa Expr
        if e.head === :call
            op = e.args[1]
            args = Any[feval(a, env) for a in e.args[2:end]]
            return fapply(op, args)
        elseif e.head === :comparison
            v = true
            for k in 1:2:length(e.args)-2
                v = v && fapply(e.args[k+1], Any[feval(e.args[k], env), feval(e.args[k+2], env)])
            end
            return v
        elseif e.head === :&& 
            return feval(e.args[1], env) && feval(e.args[2], env)
        elseif e.head === :||
            return feval(e.args[1], env) || feval(e.args[2], env)
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
    if op === :(<=) || op === :(<) || op === :(>=) || op === :(>) ||
       op === :(==) || op === :(!=)
        (args[1] isa Number && args[2] isa Number) ||
            error("unsupported: variable-dependent conditional (comparison on a model expression)")
        return getfield(Base, op)(args[1], args[2])
    elseif op === :!
        return !args[1]
    end
    f = get(JLFUN, Symbol(uppercase(String(op))), nothing)
    f === nothing && error("unsupported function $op in SIF expression")
    return f(args...)
end

"Collect call-head symbols (function names) of a parsed expression."
function callheads(e, acc = Set{Symbol}())
    if e isa Expr
        if e.head === :call && e.args[1] isa Symbol
            push!(acc, e.args[1])
        end
        for a in e.args
            callheads(a, acc)
        end
    end
    return acc
end

"A conditional assignment card in INDIVIDUALS: `target = ast` when the logical
temporary `cond` is `sense` (I-card: sense=true, E-card: sense=false)."
struct CondAssign
    cond::Symbol
    sense::Bool
    ast::Any
end

entdeps(ent) = ent isa CondAssign ? union!(freesyms(ent.ast), Set([ent.cond])) :
               freesyms(ent)
entast(ent) = ent isa CondAssign ? ent.ast : ent

"Collect free symbols of a parsed expression (call heads and comparison
operators excluded)."
function freesyms(e, acc = Set{Symbol}())
    if e isa Symbol
        (e === :true || e === :false) || push!(acc, e)
    elseif e isa Expr
        if e.head === :call
            for a in e.args[2:end]; freesyms(a, acc); end
        elseif e.head === :comparison
            for k in 1:2:length(e.args); freesyms(e.args[k], acc); end
        else
            for a in e.args
                (a isa Symbol || a isa Expr) && freesyms(a, acc)
            end
        end
    end
    return acc
end
