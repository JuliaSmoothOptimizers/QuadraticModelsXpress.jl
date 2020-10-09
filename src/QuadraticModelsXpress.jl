module QuadraticModelsXpress

export xpress

using SolverTools
using LinearAlgebra, SparseArrays
using Xpress

const xpress_statuses = Dict(0 => :unknown,
                             1 => :acceptable,
                             2 => :infeasible,
                             3 => :exception,
                             4 => :max_eval,
                             5 => :unbounded,
                             6 => :exception,
                             7 => :exception,
                             8 => :exception)

function sparse_csr(I, J, V, m=maximum(I), n=maximum(J))
    csrrowptr = zeros(Int, m+1)
    # Compute the CSR form's row counts and store them shifted forward by one in csrrowptr
    coolen = length(I)
    min(length(J), length(V)) >= coolen || throw(ArgumentError("J and V need length >= length(I) = $coolen"))
    @inbounds for k in 1:coolen
        Ik = I[k]
        if 1 > Ik || m < Ik
            throw(ArgumentError("row indices I[k] must satisfy 1 <= I[k] <= m"))
        end
        csrrowptr[Ik+1] += 1
    end

    # Compute the CSR form's rowptrs and store them shifted forward by one in csrrowptr
    countsum = 1
    csrrowptr[1] = 1
    @inbounds for i in 2:(m+1)
        overwritten = csrrowptr[i]
        csrrowptr[i] = countsum
        countsum += overwritten
    end

    # Counting-sort the column and nonzero values from J and V into csrcolval and csrnzval
    # Tracking write positions in csrrowptr corrects the row pointers
    csrcolval = zeros(Int, length(I))
    csrnzval = zeros(length(I))
    @inbounds for k in 1:coolen
        Ik, Jk = I[k], J[k]
        if 1 > Jk || n < Jk
            throw(ArgumentError("column indices J[k] must satisfy 1 <= J[k] <= n"))
        end
        csrk = csrrowptr[Ik+1]
        csrrowptr[Ik+1] = csrk + 1
        csrcolval[csrk] = Jk
        csrnzval[csrk] = V[k]
    end
    csrrowptr = csrrowptr[1:end-1]
    return csrrowptr, csrcolval, csrnzval
end


function xpress(QM; method="b", kwargs...)

    prob = Ref{Xpress.Lib.XPRSprob}()
    Xpress.Lib.XPRScreateprob(prob)
    # use kwargs change to presolve, scaling and crossover mode
    # example: xpress(QM, presolve=0) (see doc for other options)
    for (k, v) in kwargs
        if k==:presolve
            Xpress.Lib.XPRSsetintcontrol(prob.x, Xpress.Lib.XPRS_PRESOLVE, 0)  # 0 no presolve, 1=presolve
        elseif k==:scaling
            GRBsetintparam(env, "ScaleFlag", v) # 0 = no scaling
        elseif k==:crossover
            Xpress.Lib.XPRSsetintcontrol(prob.x, Xpress.Lib.XPRS_SCALING, 0)  # 0 no scaling
        elseif k==:display
            Xpress.Lib.XPRSsetintcontrol(prob.x, Xpress.Lib.XPRS_CROSSOVER, 0)  # 0 no crossover
        end
    end

    srowtypes = fill(Cchar('A'), QM.meta.ncon)
    if length(QM.meta.jinf) > 0
        error("infeasible bounds in A")
    end
    p_low, p_upp, p_rng, p_fix, p_free = 1, 1, 1, 1, 1
    for j=1:QM.meta.ncon
        if length(QM.meta.jlow) > 0 && QM.meta.jlow[p_low] == j
            srowtypes[j] = Cchar('G')
            if (p_low < length(QM.meta.jlow)) p_low += 1 end
        elseif length(QM.meta.jupp) > 0 && QM.meta.jupp[p_upp] == j
            srowtypes[j] = Cchar('L')
            if (p_upp < length(QM.meta.jupp)) p_upp += 1 end
        elseif length(QM.meta.jrng) > 0 && QM.meta.jrng[p_rng] == j
            srowtypes[j] = Cchar('R')
            if (p_rng < length(QM.meta.jrng)) p_rng += 1 end
        elseif length(QM.meta.jfix) > 0 && QM.meta.jfix[p_fix] == j
            srowtypes[j] = Cchar('E')
            if (p_fix < length(QM.meta.jfix)) p_fix += 1 end
        elseif length(QM.meta.jfree) > 0 && QM.meta.jfree[p_free] == j
            srowtypes[j] = Cchar('N')
            if (p_free < length(QM.meta.jfree)) p_free += 1 end
        else
            error("A error")
        end
    end
    rhs = zeros(QM.meta.ncon)
    drange = zeros(QM.meta.ncon)
    for j = 1:QM.meta.ncon
        if QM.meta.lcon[j] != -Inf && QM.meta.ucon[j] != Inf
            rhs[j] = QM.meta.ucon[j]
            drange[j] = QM.meta.ucon[j] - QM.meta.lcon[j]
        elseif QM.meta.lcon[j] != -Inf && QM.meta.ucon[j] == Inf
            rhs[j] = QM.meta.lcon[j]
        elseif QM.meta.lcon[j] == -Inf && QM.meta.ucon[j] != Inf
            rhs[j] = QM.meta.ucon[j]
        else
            rhs[j] = Xpress.Lib.XPRS_PLUSINFINITY
        end
    end

    A = sparse(QM.data.Arows, QM.data.Acols, QM.data.Avals, QM.meta.ncon, QM.meta.nvar)
    lvar, uvar = zeros(QM.meta.nvar), zeros(QM.meta.nvar)

    for i=1:QM.meta.nvar
        QM.meta.lvar[i] == -Inf ? lvar[i] = Xpress.Lib.XPRS_MINUSINFINITY : lvar[i] = QM.meta.lvar[i]
        QM.meta.uvar[i] == Inf ? uvar[i] = Xpress.Lib.XPRS_PLUSINFINITY : uvar[i] = QM.meta.uvar[i]
    end

    if QM.meta.nnzh > 0
        Xpress.Lib.XPRSloadqp(prob.x, QM.meta.name, QM.meta.nvar, QM.meta.ncon, srowtypes, rhs, drange,
                              QM.data.c, convert(Array{Cint,1}, A.colptr.-1), C_NULL,
                              convert(Array{Cint,1}, A.rowval.-1), A.nzval, lvar, uvar, QM.meta.nnzh,
                              convert(Array{Cint,1}, QM.data.Hrows.-1), convert(Array{Cint,1}, QM.data.Hcols.-1),
                              QM.data.Hvals)
    else
        Xpress.Lib.XPRSloadlp(prob.x, QM.meta.name, QM.meta.nvar, QM.meta.ncon, srowtypes, rhs, drange,
                              QM.data.c, convert(Array{Cint,1}, A.colptr.-1), C_NULL,
                              convert(Array{Cint,1}, A.rowval.-1), A.nzval, lvar, uvar)
    end
    Xpress.Lib.XPRSchgobj(prob.x, Cint(1), Cint.([-1]), [-QM.data.c0])

    t = @timed begin
        Xpress.Lib.XPRSlpoptimize(prob.x, method)
    end

    x, y, s = zeros(QM.meta.nvar), zeros(QM.meta.ncon), zeros(QM.meta.ncon), zeros(QM.meta.nvar)
    Xpress.Lib.XPRSgetsol(prob.x, x, C_NULL, y, s)
    baritcnt = Ref{Cint}()
    Xpress.Lib.XPRSgetintattrib(prob.x, Xpress.Lib.XPRS_BARITER, baritcnt)
    objval = Ref{Float64}()
    Xpress.Lib.XPRSgetdblattrib(prob.x, Xpress.Lib.XPRS_LPOBJVAL, objval)
    status = Ref{Cint}()
    Xpress.Lib.XPRSgetintattrib(prob.x, Xpress.Lib.XPRS_LPSTATUS, status)
    p_feas = Ref{Float64}() # sum
    Xpress.Lib.XPRSgetdblattrib(prob.x, Xpress.Lib.XPRS_BARPRIMALINF, p_feas)
    d_feas = Ref{Float64}() # sum
    Xpress.Lib.XPRSgetdblattrib(prob.x, Xpress.Lib.XPRS_BARDUALINF, d_feas)
    stats = GenericExecutionStats(get(xpress_statuses, status[], :unknown),
                                  QM, solution = x,
                                  objective = objval[],
                                  primal_feas = p_feas[],
                                  dual_feas = d_feas[],
                                  iter = Int64(baritcnt[]),
                                  multipliers = y,
                                  elapsed_time = t[2])
    return stats
end

end
