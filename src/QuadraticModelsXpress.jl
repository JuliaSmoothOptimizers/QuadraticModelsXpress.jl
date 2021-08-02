module QuadraticModelsXpress

export xpress

using QuadraticModels, SolverCore
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

function xpress(QM::QuadraticModel; method="b", kwargs...)

    prob = Xpress.XpressProblem()
    # use kwargs change to presolve, scaling and crossover mode
    # example: xpress(QM, presolve=0, bargapstop=1e-8) 
    for (k, v) in kwargs
        if k==:presolve
            Xpress.setintcontrol(prob, Xpress.Lib.XPRS_PRESOLVE, v)  # 0 no presolve, 1=presolve
        elseif k==:scaling
            Xpress.setintcontrol(prob, Xpress.Lib.XPRS_SCALING, v)  # 0 no scaling
        elseif k==:crossover
            Xpress.setintcontrol(prob, Xpress.Lib.XPRS_CROSSOVER, v)  # 0 no crossover
        elseif k==:bargapstop
            Xpress.setdblcontrol(prob, Xpress.Lib.XPRS_BARGAPSTOP, v)
        elseif k==:barprimalstop
            Xpress.setdblcontrol(prob, Xpress.Lib.XPRS_BARPRIMALSTOP, v)
        elseif k==:bardualstop
            Xpress.setdblcontrol(prob, Xpress.Lib.XPRS_BARDUALSTOP, v)
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
        Xpress.loadqp(prob, QM.meta.name, QM.meta.nvar, QM.meta.ncon, srowtypes, rhs, drange,
                      QM.data.c, convert(Array{Cint,1}, A.colptr.-1), C_NULL,
                      convert(Array{Cint,1}, A.rowval.-1), A.nzval, lvar, uvar, QM.meta.nnzh,
                      convert(Array{Cint,1}, QM.data.Hrows.-1), convert(Array{Cint,1}, QM.data.Hcols.-1),
                      QM.data.Hvals)
    else
        Xpress.loadlp(prob, "", QM.meta.nvar, QM.meta.ncon, srowtypes, rhs, drange,
                      QM.data.c,convert(Array{Cint,1}, A.colptr.-1), C_NULL,convert(Array{Cint,1}, A.rowval.-1),
                      A.nzval, lvar, uvar)
    end
    Xpress.chgobj(prob, [0], [-QM.data.c0])

    start_time = time()
    Xpress.lpoptimize(prob, method)
    elapsed_time = time() - start_time

    x, y, s = zeros(QM.meta.nvar), zeros(QM.meta.ncon), zeros(QM.meta.nvar)
    Xpress.getsol(prob, x, C_NULL, y, s)
    baritcnt = Xpress.getintattrib(prob, Xpress.Lib.XPRS_BARITER)
    objval = Xpress.getdblattrib(prob, Xpress.Lib.XPRS_LPOBJVAL)
    status = Xpress.getintattrib(prob, Xpress.Lib.XPRS_LPSTATUS)
    p_feas = Xpress.getdblattrib(prob, Xpress.Lib.XPRS_BARPRIMALINF)
    d_feas = Xpress.getdblattrib(prob, Xpress.Lib.XPRS_BARDUALINF)

    Xpress.destroyprob(prob)

    stats = GenericExecutionStats(get(xpress_statuses, status, :unknown),
                                  QM, solution = x,
                                  objective = objval,
                                  primal_feas = p_feas,
                                  dual_feas = d_feas,
                                  iter = Int64(baritcnt),
                                  multipliers = y,
                                  elapsed_time = elapsed_time)
    return stats
end

end
