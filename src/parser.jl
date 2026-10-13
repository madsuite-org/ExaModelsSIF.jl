# parser.jl — SIF file front-end: fixed-field card tokenizer + data-section
# interpreter (DO loops, parameter arithmetic, array name expansion) producing
# a fully instantiated problem IR, plus the ELEMENTS/GROUPS function-part
# parser. Clean-room implementation from the published SIF reference
# (https://www.numerical.rl.ac.uk/lancelot/sif/); no CUTEst code consulted.

# ------------------------------------------------------------------ IR types

struct ElType
    name::String
    evars::Vector{Symbol}                 # elemental variables, declaration order
    ivars::Vector{Symbol}                 # internal variables (empty if none)
    epars::Vector{Symbol}                 # element parameters, declaration order
    rrows::Dict{Symbol,Vector{Pair{Symbol,Float64}}}  # ivar => [(evar => coef)]
    temps::Vector{Pair{Symbol,Any}}       # ordered A-assignments (name => ast)
    fexpr::Any                            # F expression ast
end

struct GType
    name::String
    gvar::Symbol
    gpars::Vector{Symbol}
    temps::Vector{Pair{Symbol,Any}}
    fexpr::Any
end

mutable struct SIFElement
    name::String
    etype::String                          # "" until typed
    varmap::Dict{Symbol,Int}               # elemental var -> problem var index
    parvals::Dict{Symbol,Float64}
end

mutable struct SIFGroup
    name::String
    kind::Char                             # 'N' objective, 'E','L','G' constraint
    gtype::Union{Nothing,String}           # nothing => use default/trivial
    scale::Float64
    constant::Float64
    range_::Union{Nothing,Float64}
    lin::Vector{Pair{Int,Float64}}         # var index => coef (duplicates summed later)
    elems::Vector{Tuple{String,Float64}}   # (element name, weight)
    gparvals::Dict{Symbol,Float64}
end

mutable struct SIFProblem
    name::String
    varnames::Vector{String}
    varidx::Dict{String,Int}
    lvar::Vector{Float64}
    uvar::Vector{Float64}
    x0::Vector{Float64}
    varscale::Vector{Float64}
    grouporder::Vector{String}
    groups::Dict{String,SIFGroup}
    elemorder::Vector{String}
    elements::Dict{String,SIFElement}
    eltypes::Dict{String,ElType}
    gtypes::Dict{String,GType}
    default_etype::Union{Nothing,String}
    default_gtype::Union{Nothing,String}
    quad::Vector{Tuple{Int,Int,Float64}}   # (i, j, coef)
end

SIFProblem() = SIFProblem("", String[], Dict{String,Int}(),
    Float64[], Float64[], Float64[], Float64[],
    String[], Dict{String,SIFGroup}(),
    String[], Dict{String,SIFElement}(),
    Dict{String,ElType}(), Dict{String,GType}(),
    nothing, nothing, Tuple{Int,Int,Float64}[])

# ------------------------------------------------------------------ tokenizer

struct Card
    code::String
    f2::String
    f3::String
    f4::String
    f5::String
    f6::String
    expr::String   # raw text from col 25 to EOL (for A/F/G/H lines)
end

function card(ln::AbstractString)
    di = findfirst(" \$", ln)
    di !== nothing && (ln = ln[1:first(di)-1])
    L = lastindex(ln)
    g(a, b) = a <= L ? String(strip(ln[a:min(b, L)])) : ""
    # (name, value) pairs: the name is the 10-column name subfield verbatim
    # (names may contain interior blanks, e.g. HELSBY "N      1") and the value
    # is the rest — EXCEPT when the name field's trailing token joined with the
    # value field forms a number (misaligned value bleeding left, TORSIONA
    # "C(I,J)   0" + ".25" -> name "C(I,J)", value "0.25").
    function pair(a, b)
        a > L && return ("", "")
        nfraw = ln[a:min(a + 9, L)]
        vf = a + 10 <= L ? String(strip(ln[a+10:min(b, L)])) : ""
        # Merge a value bleeding across the subfield boundary (TORSION* class)
        # ONLY when the boundary itself is contiguous non-blank text; names may
        # legitimately contain interior blanks (HELSBY "N      1").
        if length(nfraw) == 10 && !isspace(nfraw[end]) &&
           a + 10 <= L && !isspace(ln[a+10])
            mm = match(r"^(.*\S)\s+(\S+)$", nfraw)
            if mm !== nothing
                merged = replace(mm[2] * vf, " " => "")
                if tryparse(Float64, replace(merged, r"[dD](?=[0-9+-])" => "e")) !== nothing
                    return (String(mm[1]), merged)
                end
            end
        end
        return (String(strip(nfraw)), vf)
    end
    f3, f4 = pair(15, 36)
    f5, f6 = pair(40, 61)
    Card(g(2, 3), g(5, 14), f3, f4, f5, f6,
         L >= 25 ? String(rstrip(ln[25:L])) : "")
end

issectionline(ln) = !isempty(ln) && !isspace(ln[1]) && ln[1] != '*'

# Split raw lines into the data part and the function parts.
function splitparts(lines)
    datacards = Any[]                       # Card or (:section, name)
    elemlines = String[]; grouplines = String[]
    part = :data
    for ln in lines
        isempty(strip(ln)) && continue
        startswith(ln, '*') && continue
        if issectionline(ln)
            toks = split(ln)
            head = uppercase(String(toks[1]))
            if part === :data
                if head == "ENDATA"
                    part = :pending
                elseif head == "NAME"
                    push!(datacards, (:name, length(toks) > 1 ? String(toks[2]) : ""))
                else
                    sec = join(toks, " ")  # e.g. "START POINT", may include trailing junk
                    push!(datacards, (:section, uppercase(sec)))
                end
            else
                if head == "ELEMENTS"
                    part = :elements
                elseif head == "GROUPS"
                    part = :groups
                elseif head == "ENDATA"
                    part = :pending
                elseif part === :elements
                    push!(elemlines, ln)
                elseif part === :groups
                    push!(grouplines, ln)
                end
            end
        else
            if part === :data
                push!(datacards, card(ln))
            elseif part === :elements
                push!(elemlines, ln)
            elseif part === :groups
                push!(grouplines, ln)
            end
        end
    end
    return datacards, elemlines, grouplines
end

# ------------------------------------------------------------- interpreter

mutable struct State
    P::SIFProblem
    section::String
    ip::Dict{String,Int}
    rp::Dict{String,Float64}
    loopfinal::Symbol          # :last or :fortran
    firstset::Dict{String,String}  # section => first set name seen (field 2)
end

parsenum(s::AbstractString) =
    parse(Float64, replace(s, r"[dD](?=[0-9+-])" => "e", " " => ""))

function resolveint(st::State, tok::AbstractString)
    v = tryparse(Int, tok)
    v !== nothing && return v
    haskey(st.ip, tok) && return st.ip[tok]
    if occursin('(', tok)
        t = expandname(st, tok)
        haskey(st.ip, t) && return st.ip[t]
        v = tryparse(Int, t)
        v !== nothing && return v
    end
    error("undefined integer parameter '$tok'")
end

"Expand an array name like X(I+1) or Y(1,J) to its flat instantiated name."
function expandname(st::State, name::AbstractString)
    occursin('(', name) || return String(name)
    return replace(name, r"\(([^()]*)\)" => m -> begin
        inner = m[2:end-1]
        vals = String[]
        for tok in split(inner, ',')
            v = try
                resolveint(st, strip(tok))
            catch
                nothing
            end
            v === nothing && return m   # not an index list: part of the name
            push!(vals, string(v))
        end
        join(vals, ",")
    end)
end

getrp(st::State, tok::AbstractString) = begin
    t = expandname(st, tok)
    haskey(st.rp, t) || error("undefined real parameter '$t'")
    st.rp[t]
end

const IPARAM = Set(["IE", "IA", "IS", "IM", "ID", "I=", "I+", "I-", "I*", "I/", "IR"])
const RPARAM = Set(["RE", "RA", "RS", "RM", "RD", "R=", "R+", "R-", "R*", "R/",
                    "RI", "RF", "R("])
# Array real-parameter codes: identical semantics to the R-codes, with array
# names (expanded by expandname like any other name).
const APARAM = Dict("AE" => "RE", "AA" => "RA", "AS" => "RS", "AM" => "RM",
    "AD" => "RD", "A=" => "R=", "A+" => "R+", "A-" => "R-", "A*" => "R*",
    "A/" => "R/", "AI" => "RI", "AF" => "RF", "A(" => "R(")

function handle_iparam!(st::State, cd::Card)
    dest = expandname(st, cd.f2)
    c = cd.code
    if c in ("IA", "IS", "IM", "ID") && isempty(cd.f4)
        c = c == "IA" ? "I+" : c == "IS" ? "I-" : c == "IM" ? "I*" : "I/"
    end
    if c == "IE"
        st.ip[dest] = resolveint(st, isempty(cd.f4) ? cd.f3 : cd.f4)
    elseif c == "IA"
        st.ip[dest] = st.ip[expandname(st, cd.f3)] + resolveint(st, cd.f4)
    elseif c == "IS"
        st.ip[dest] = resolveint(st, cd.f4) - st.ip[expandname(st, cd.f3)]
    elseif c == "IM"
        st.ip[dest] = st.ip[expandname(st, cd.f3)] * resolveint(st, cd.f4)
    elseif c == "ID"
        st.ip[dest] = div(resolveint(st, cd.f4), st.ip[expandname(st, cd.f3)])
    elseif c == "I="
        st.ip[dest] = st.ip[expandname(st, cd.f3)]
    elseif c == "I+"
        st.ip[dest] = st.ip[expandname(st, cd.f3)] + st.ip[expandname(st, cd.f5)]
    elseif c == "I-"
        st.ip[dest] = st.ip[expandname(st, cd.f3)] - st.ip[expandname(st, cd.f5)]
    elseif c == "I*"
        st.ip[dest] = st.ip[expandname(st, cd.f3)] * st.ip[expandname(st, cd.f5)]
    elseif c == "I/"
        st.ip[dest] = div(st.ip[expandname(st, cd.f3)], st.ip[expandname(st, cd.f5)])
    elseif c == "IR"
        st.ip[dest] = trunc(Int, st.rp[expandname(st, cd.f3)])
    end
end

function handle_rparam!(st::State, cd::Card)
    dest = expandname(st, cd.f2)
    c = cd.code
    # RA/RS/RM/RD with an empty value field take a PARAMETER operand in field 5
    # and act like their two-parameter cousins (f3 op f5), e.g.
    # "RD DP/NI  DP  RNI" means DP/NI := DP / RNI.
    if c in ("RA", "RS", "RM", "RD") && isempty(cd.f4)
        c = c == "RA" ? "R+" : c == "RS" ? "R-" : c == "RM" ? "R*" : "R/"
    end
    if c == "RE"
        st.rp[dest] = parsenum(isempty(cd.f4) ? cd.f3 : cd.f4)
    elseif c == "RA"
        st.rp[dest] = st.rp[expandname(st, cd.f3)] + parsenum(cd.f4)
    elseif c == "RS"
        st.rp[dest] = parsenum(cd.f4) - st.rp[expandname(st, cd.f3)]
    elseif c == "RM"
        st.rp[dest] = st.rp[expandname(st, cd.f3)] * parsenum(cd.f4)
    elseif c == "RD"
        st.rp[dest] = parsenum(cd.f4) / st.rp[expandname(st, cd.f3)]
    elseif c == "R="
        st.rp[dest] = st.rp[expandname(st, cd.f3)]
    elseif c == "R+"
        st.rp[dest] = st.rp[expandname(st, cd.f3)] + st.rp[expandname(st, cd.f5)]
    elseif c == "R-"
        st.rp[dest] = st.rp[expandname(st, cd.f3)] - st.rp[expandname(st, cd.f5)]
    elseif c == "R*"
        st.rp[dest] = st.rp[expandname(st, cd.f3)] * st.rp[expandname(st, cd.f5)]
    elseif c == "R/"
        st.rp[dest] = st.rp[expandname(st, cd.f3)] / st.rp[expandname(st, cd.f5)]
    elseif c == "RI"
        st.rp[dest] = Float64(st.ip[expandname(st, cd.f3)])
    elseif c == "RF"
        st.rp[dest] = rfun(cd.f3)(parsenum(cd.f4))
    elseif c == "R("
        st.rp[dest] = rfun(cd.f3)(getrp(st, cd.f5))
    end
end

function rfun(name)
    f = get(JLFUN, Symbol(uppercase(String(name))), nothing)
    f === nothing && error("unsupported parameter function '$name'")
    return f
end

# -------- problem-object helpers

function getvar!(st::State, name::AbstractString)
    P = st.P
    haskey(P.varidx, name) && return P.varidx[name]
    push!(P.varnames, String(name))
    push!(P.lvar, 0.0)        # SIF default bounds: [0, +Inf)
    push!(P.uvar, Inf)
    push!(P.x0, 0.0)
    push!(P.varscale, 1.0)
    P.varidx[name] = length(P.varnames)
end

function getgroup!(st::State, name::AbstractString, kind::Union{Nothing,Char})
    P = st.P
    if haskey(P.groups, name)
        g = P.groups[name]
        if kind !== nothing && (g.kind == 'N') != (kind == 'N')
            alias = name * "\x01old"
            haskey(P.groups, alias) && error("group $name redeclared twice across namespaces")
            P.groups[alias] = g
            i = findfirst(==(name), P.grouporder)
            P.grouporder[i] = alias
            g.name = alias
        else
            kind !== nothing && g.kind != kind &&
                error("group $name redeclared with kind $kind (was $(g.kind))")
            return g
        end
    end
    kind === nothing && error("group '$name' referenced before declaration")
    g = SIFGroup(String(name), kind, nothing, 1.0, 0.0, nothing,
                 Pair{Int,Float64}[], Tuple{String,Float64}[], Dict{Symbol,Float64}())
    P.groups[name] = g
    push!(P.grouporder, String(name))
    return g
end

function getelem!(st::State, name::AbstractString)
    P = st.P
    haskey(P.elements, name) && return P.elements[name]
    e = SIFElement(String(name), "", Dict{Symbol,Int}(), Dict{Symbol,Float64}())
    P.elements[name] = e
    push!(P.elemorder, String(name))
    return e
end

unquote(s) = startswith(s, '\'') && endswith(s, '\'') ? s[2:end-1] : s

# -------- section card handlers

function handle_card!(st::State, cd::Card)
    code = cd.code
    code in IPARAM && return handle_iparam!(st, cd)
    code in RPARAM && return handle_rparam!(st, cd)
    haskey(APARAM, code) && return handle_rparam!(st,
        Card(APARAM[code], cd.f2, cd.f3, cd.f4, cd.f5, cd.f6, cd.expr))
    sec = st.section
    if sec == "" && isempty(code)
        return   # stray unmarked comment line before the first section
    end
    # Only the first named set applies in these sections; later sets are
    # alternative constants/bounds/starts (e.g. recorded solutions).
    if sec in ("CONSTANTS", "RHS", "RHS'", "RANGES", "BOUNDS", "START POINT") &&
       !isempty(cd.f2)
        fs = get!(st.firstset, sec, cd.f2)
        fs == cd.f2 || return
    end
    if sec == "VARIABLES" || sec == "COLUMNS"
        handle_variables!(st, cd)
    elseif sec == "GROUPS" || sec == "ROWS" || sec == "CONSTRAINTS"
        handle_groups!(st, cd)
    elseif sec == "CONSTANTS" || sec == "RHS" || sec == "RHS'"
        handle_constants!(st, cd)
    elseif sec == "RANGES"
        handle_ranges!(st, cd)
    elseif sec == "BOUNDS"
        handle_bounds!(st, cd)
    elseif sec == "START POINT"
        handle_start!(st, cd)
    elseif sec in ("QUADRATIC", "QUADS", "QUADOBJ", "HESSIAN", "QMATRIX", "QSECTION")
        handle_quadratic!(st, cd)
    elseif sec == "ELEMENT TYPE"
        handle_eltype!(st, cd)
    elseif sec == "ELEMENT USES"
        handle_eluses!(st, cd)
    elseif sec == "GROUP TYPE"
        handle_gtype!(st, cd)
    elseif sec == "GROUP USES"
        handle_guses!(st, cd)
    elseif sec == "OBJECT BOUND"
        return
    else
        error("card in unsupported section '$sec' (code '$code')")
    end
end

function handle_variables!(st::State, cd::Card)
    cd.code in ("", "X", "XN", "Z") ||
        error("unsupported VARIABLES code '$(cd.code)'")
    v = getvar!(st, expandname(st, cd.f2))
    if unquote(cd.f3) == "SCALE"
        # Variable scale factors are LANCELOT hints; SIFdecode/CUTEst evaluate
        # the UNSCALED problem (calibrated on MEYER3: ref x0/obj are the raw
        # SIF values). Parsed and discarded.
        return
    end
    if cd.code == "Z"
        isempty(cd.f3) && return
        g = getgroup!(st, expandname(st, cd.f3), nothing)
        push!(g.lin, v => getrp(st, cd.f5))
        return
    end
    # COLUMNS-style (group, coef) pairs; INTEGER/CONTINUOUS markers are
    # accepted and ignored (continuous relaxation, as CUTEst evaluates it)
    for (gn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
        isempty(gn) && continue
        unquote(gn) in ("INTEGER", "CONTINUOUS") && continue
        g = getgroup!(st, expandname(st, gn), nothing)
        push!(g.lin, v => isempty(val) ? 1.0 : parsenum(val))
    end
end

function handle_groups!(st::State, cd::Card)
    code = cd.code
    kindc = code[end]
    kindc in ('N', 'E', 'G', 'L') ||
        error("unsupported GROUPS code '$code'")
    prefix = length(code) == 2 ? code[1] : ' '
    g = getgroup!(st, expandname(st, cd.f2), kindc)
    if unquote(cd.f3) == "SCALE"
        g.scale = prefix == 'Z' ? getrp(st, cd.f5) : parsenum(cd.f4)
        return
    end
    if prefix == 'Z'
        isempty(cd.f3) && return
        v = getvar!(st, expandname(st, cd.f3))
        push!(g.lin, v => getrp(st, cd.f5))
    else
        for (vn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
            isempty(vn) && continue
            v = getvar!(st, expandname(st, vn))
            push!(g.lin, v => isempty(val) ? 1.0 : parsenum(val))
        end
    end
end

function handle_constants!(st::State, cd::Card)
    code = cd.code
    if startswith(code, "Z")
        g = getgroup!(st, expandname(st, cd.f3), nothing)
        g.constant = getrp(st, cd.f5)
        return
    end
    for (gn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
        isempty(gn) && continue
        if unquote(gn) == "DEFAULT"
            v = parsenum(val)
            for g in values(st.P.groups); g.constant = v; end
        else
            g = getgroup!(st, expandname(st, gn), nothing)
            g.constant = parsenum(val)
        end
    end
end

function handle_ranges!(st::State, cd::Card)
    code = cd.code
    if startswith(code, "Z")
        g = getgroup!(st, expandname(st, cd.f3), nothing)
        g.range_ = getrp(st, cd.f5)
        return
    end
    for (gn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
        isempty(gn) && continue
        if unquote(gn) == "DEFAULT"
            v = parsenum(val)
            for g in values(st.P.groups); g.range_ = v; end
        else
            g = getgroup!(st, expandname(st, gn), nothing)
            g.range_ = parsenum(val)
        end
    end
end

function eachboundvar(st::State, tok)
    P = st.P
    return unquote(tok) == "DEFAULT" ? (1:length(P.varnames)) :
        (P.varidx[expandname(st, tok)]:P.varidx[expandname(st, tok)])
end

function handle_bounds!(st::State, cd::Card)
    code = cd.code
    P = st.P
    vs = eachboundvar(st, cd.f3)
    if code in ("FR", "XR", "ZR")
        for v in vs; P.lvar[v] = -Inf; P.uvar[v] = Inf; end
    elseif code in ("MI", "XM", "ZM")
        for v in vs; P.lvar[v] = -Inf; end
    elseif code in ("PL", "XP", "ZP")
        for v in vs; P.uvar[v] = Inf; end
    elseif code in ("LO", "XL", "ZL", "UP", "XU", "ZU", "FX", "XX", "ZX")
        val = startswith(code, "Z") ? getrp(st, cd.f5) : parsenum(cd.f4)
        if code in ("LO", "XL", "ZL")
            for v in vs; P.lvar[v] = val; end
        elseif code in ("UP", "XU", "ZU")
            for v in vs; P.uvar[v] = val; end
        else
            for v in vs; P.lvar[v] = val; P.uvar[v] = val; end
        end
    else
        error("unsupported BOUNDS code '$code'")
    end
end

function handle_start!(st::State, cd::Card)
    code = cd.code
    P = st.P
    code in ("", "X", "V", "XV", "ZV", "Z") ||
        return   # dual-value starts and other exotics: ignore
    if code == "ZV" || code == "Z"
        tok = cd.f3
        val = getrp(st, cd.f5)
        if unquote(tok) == "DEFAULT"
            fill!(P.x0, val)
        else
            P.x0[P.varidx[expandname(st, tok)]] = val
        end
        return
    end
    for (tok, sval) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
        isempty(tok) && continue
        val = isempty(sval) ? 0.0 : parsenum(sval)
        if unquote(tok) == "DEFAULT"
            fill!(P.x0, val)
        else
            P.x0[P.varidx[expandname(st, tok)]] = val
        end
    end
end

function handle_quadratic!(st::State, cd::Card)
    code = cd.code
    i = getvar!(st, expandname(st, cd.f2))
    if startswith(code, "Z")
        j = getvar!(st, expandname(st, cd.f3))
        push!(st.P.quad, (i, j, getrp(st, cd.f5)))
        return
    end
    for (vn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
        isempty(vn) && continue
        j = getvar!(st, expandname(st, vn))
        push!(st.P.quad, (i, j, parsenum(val)))
    end
end

# element/group type DECLARATIONS are accumulated here, their function bodies
# attached later by the function-part parser.
mutable struct TypeDecl
    evars::Vector{Symbol}
    ivars::Vector{Symbol}
    epars::Vector{Symbol}
end
const _typedecls = Dict{String,TypeDecl}()   # reset per parse (stored in State? kept simple)

function handle_eltype!(st::State, cd::Card)
    code = cd.code
    d = get!(_typedecls, cd.f2, TypeDecl(Symbol[], Symbol[], Symbol[]))
    tgt = code == "EV" ? d.evars : code == "IV" ? d.ivars : code == "EP" ? d.epars :
        error("unsupported ELEMENT TYPE code '$code'")
    for tok in (cd.f3, cd.f5)
        isempty(tok) && continue
        push!(tgt, Symbol(tok))
    end
end

function handle_eluses!(st::State, cd::Card)
    code = cd.code
    if code in ("T", "XT")
        if unquote(cd.f2) == "DEFAULT"
            st.P.default_etype = cd.f3
        else
            e = getelem!(st, expandname(st, cd.f2))
            e.etype = cd.f3
        end
    elseif code in ("V", "XV", "ZV")
        e = getelem!(st, expandname(st, cd.f2))
        e.varmap[Symbol(cd.f3)] = getvar!(st, expandname(st, cd.f5))
    elseif code in ("P", "XP")
        e = getelem!(st, expandname(st, cd.f2))
        for (pn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
            isempty(pn) && continue
            e.parvals[Symbol(pn)] = parsenum(val)
        end
    elseif code == "ZP"
        e = getelem!(st, expandname(st, cd.f2))
        e.parvals[Symbol(cd.f3)] = getrp(st, cd.f5)
    else
        error("unsupported ELEMENT USES code '$code'")
    end
end

function handle_gtype!(st::State, cd::Card)
    code = cd.code
    d = get!(_typedecls, "G#" * cd.f2, TypeDecl(Symbol[], Symbol[], Symbol[]))
    if code == "GV"
        push!(d.evars, Symbol(cd.f3))
    elseif code == "GP"
        for tok in (cd.f3, cd.f5)
            isempty(tok) && continue
            push!(d.epars, Symbol(tok))
        end
    else
        error("unsupported GROUP TYPE code '$code'")
    end
end

function handle_guses!(st::State, cd::Card)
    code = cd.code
    if isempty(code) && haskey(_typedecls, "G#" * cd.f3)
        code = "T"
    end
    if code in ("T", "XT")
        if unquote(cd.f2) == "DEFAULT"
            st.P.default_gtype = cd.f3
        else
            g = getgroup!(st, expandname(st, cd.f2), nothing)
            g.gtype = cd.f3
        end
    elseif code in ("E", "XE")
        g = getgroup!(st, expandname(st, cd.f2), nothing)
        for (en, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
            isempty(en) && continue
            w = isempty(val) ? 1.0 : parsenum(val)
            push!(g.elems, (expandname(st, en), w))
        end
    elseif code == "ZE"
        g = getgroup!(st, expandname(st, cd.f2), nothing)
        push!(g.elems, (expandname(st, cd.f3), getrp(st, cd.f5)))
    elseif code in ("P", "XP")
        g = getgroup!(st, expandname(st, cd.f2), nothing)
        for (pn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
            isempty(pn) && continue
            g.gparvals[Symbol(pn)] = parsenum(val)
        end
    elseif code == "ZP"
        g = getgroup!(st, expandname(st, cd.f2), nothing)
        g.gparvals[Symbol(cd.f3)] = getrp(st, cd.f5)
    else
        error("unsupported GROUP USES code '$code'")
    end
end

# -------- card list execution with DO loops

function run_cards!(st::State, cards, i0::Int, i1::Int)
    i = i0
    while i <= i1
        cd = cards[i]
        if cd isa Tuple
            if cd[1] === :section
                st.section = cd[2]
            elseif cd[1] === :name
                st.P.name = cd[2]
            end
            i += 1
            continue
        end
        if cd.code == "DO"
            var = cd.f2
            lo = resolveint(st, cd.f3)
            hi = resolveint(st, cd.f5)
            step = 1
            j = i + 1
            while j <= i1 && cards[j] isa Card && cards[j].code == "DI"
                step = resolveint(st, cards[j].f3)
                j += 1
            end
            # OD ends the innermost loop; ND ends ALL open loops, so it
            # terminates this scan at any depth. An inner loop whose shared ND
            # lies outside the enclosing body slice terminates at end-of-range.
            depth = 1
            k = j
            while k <= i1
                c2 = cards[k]
                if c2 isa Card
                    c2.code == "DO" && (depth += 1)
                    if c2.code == "OD"
                        depth -= 1
                        depth == 0 && break
                    elseif c2.code == "ND"
                        break
                    end
                end
                k += 1
            end
            v = lo
            while (step > 0 ? v <= hi : v >= hi)
                st.ip[var] = v
                run_cards!(st, cards, j, k - 1)
                v += step
            end
            # loop variable value after completion (calibrated vs SIFDecode)
            st.ip[var] = st.loopfinal === :fortran ? v : hi
            i = k + 1
        elseif cd.code == "ND" || cd.code == "OD" || cd.code == "DI"
            error("loop control card outside loop")
        else
            try
                handle_card!(st, cd)
            catch err
                err isa ErrorException && startswith(err.msg, "[card") && rethrow()
                error("[card '$(cd.code)' f2='$(cd.f2)' f3='$(cd.f3)' f4='$(cd.f4)' " *
                      "f5='$(cd.f5)' f6='$(cd.f6)' sec='$(st.section)'] " *
                      sprint(showerror, err))
            end
            i += 1
        end
    end
end

# -------- function-part (ELEMENTS / GROUPS) parser

function parse_funpart(lines, isgroup::Bool)
    fdecls = Set{Symbol}()
    temps_global = Pair{Symbol,Any}[]
    blocks = Dict{String,NamedTuple}()
    cur = nothing            # current type name
    curtemps = Pair{Symbol,Any}[]
    currows = Dict{Symbol,Vector{Pair{Symbol,Float64}}}()
    curf = nothing
    section = ""
    pend_code = ""           # for continuation lines
    pend_name = Symbol("")
    pend_expr = ""
    pend2 = Symbol("")       # condition temp for I/E cards

    function flushpend!()
        isempty(pend_code) && return
        if pend_code == "A"
            push!(section == "GLOBALS" ? temps_global : curtemps,
                  pend_name => fparse(pend_expr))
        elseif pend_code == "F"
            curf = fparse(pend_expr)
        elseif pend_code == "I" || pend_code == "E"
            push!(curtemps, pend_name =>
                CondAssign(pend2, pend_code == "I", fparse(pend_expr)))
        end
        pend_code = ""
    end
    function flushblock!()
        flushpend!()
        cur === nothing && return
        blocks[cur] = (temps = curtemps, rows = currows, fexpr = curf)
        cur = nothing
    end

    for ln in lines
        if issectionline(ln)
            flushpend!()
            section = uppercase(String(split(ln)[1]))
            continue
        end
        cd = card(ln)
        code = cd.code
        if endswith(code, "+") && !isempty(pend_code)
            pend_expr *= " " * cd.expr
            continue
        end
        flushpend!()
        if section == "TEMPORARIES"
            code == "F" && push!(fdecls, Symbol(cd.f2))
            continue  # R/I/M/L declarations carry no information we need
        elseif section == "GLOBALS"
            if code == "A"
                pend_code = "A"; pend_name = Symbol(cd.f2); pend_expr = cd.expr
            end
        elseif section == "INDIVIDUALS"
            if code == "T"
                flushblock!()
                cur = cd.f2
                curtemps = Pair{Symbol,Any}[]
                currows = Dict{Symbol,Vector{Pair{Symbol,Float64}}}()
                curf = nothing
            elseif code == "A"
                pend_code = "A"; pend_name = Symbol(cd.f2); pend_expr = cd.expr
            elseif code == "F"
                pend_code = "F"; pend_expr = cd.expr
            elseif code == "G" || code == "H"
                pend_code = "G"; pend_expr = cd.expr   # parsed then discarded
            elseif code == "I" || code == "E"
                pend_code = code
                pend2 = Symbol(cd.f2)
                pend_name = Symbol(cd.f3)
                # expression begins after the two name fields (col 25+); strip
                # the condition/target tokens that the expr span includes
                pend_expr = cd.expr
            elseif code == "R"
                rows = get!(currows, Symbol(cd.f2), Pair{Symbol,Float64}[])
                for (vn, val) in ((cd.f3, cd.f4), (cd.f5, cd.f6))
                    isempty(vn) && continue
                    push!(rows, Symbol(vn) => parsenum(val))
                end
            end
        end
    end
    flushblock!()
    return temps_global, blocks, fdecls
end

# ------------------------------------------------------------------ driver

"""
    sif_decode(path; loopfinal = :last) -> SIFProblem

Parse and interpret a SIF file into a fully instantiated `SIFProblem` IR.
"""
function sif_decode(path::AbstractString; loopfinal::Symbol = :last)
    lines = readlines(path)
    datacards, elemlines, grouplines = splitparts(lines)
    empty!(_typedecls)
    st = State(SIFProblem(), "", Dict{String,Int}(), Dict{String,Float64}(),
               loopfinal, Dict{String,String}())
    run_cards!(st, datacards, 1, length(datacards))
    P = st.P

    gtemps_e, eblocks, fdecl_e = parse_funpart(elemlines, false)
    gtemps_g, gblocks, fdecl_g = parse_funpart(grouplines, true)

    function checkexternal(tname, temps, fexpr, fdecls)
        isempty(fdecls) && return
        syms = fexpr === nothing ? Set{Symbol}() : union!(freesyms(fexpr), callheads(fexpr))
        for (nm, ent) in temps
            union!(syms, entdeps(ent))
            union!(syms, callheads(entast(ent)))
        end
        for fn in fdecls
            fn in syms && throw(ExternalFortranError(String(tname), fn))
        end
    end

    for (tname, blk) in eblocks
        d = get(_typedecls, tname, nothing)
        d === nothing && error("element type $tname has a function body but no declaration")
        temps = vcat(gtemps_e, blk.temps)
        checkexternal(tname, temps, blk.fexpr, fdecl_e)
        P.eltypes[tname] = ElType(tname, d.evars, d.ivars, d.epars, blk.rows,
                                  temps, blk.fexpr)
    end
    for (tname, blk) in gblocks
        d = get(_typedecls, "G#" * tname, nothing)
        d === nothing && error("group type $tname has a function body but no declaration")
        isempty(d.evars) && error("group type $tname has no GV variable")
        temps = vcat(gtemps_g, blk.temps)
        checkexternal(tname, temps, blk.fexpr, fdecl_g)
        P.gtypes[tname] = GType(tname, d.evars[1], d.epars,
                                temps, blk.fexpr)
    end
    return P
end
