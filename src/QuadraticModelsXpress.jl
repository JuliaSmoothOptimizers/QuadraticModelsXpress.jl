module QuadraticModelsXpress

export xpress

using QuadraticModels
using SparseArrays
using Xpress


function xpress(QM)
    Xpress.init()
    prob = Xpress.XpressProblem()
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

    Xpress.loadlp(prob, QM.meta.name, QM.meta.nvar, QM.meta.ncon, srowtypes,
         rhs, drange, QM.data.c, A.colptr.-1, C_NULL, A.rowval.-1, A.nzval, lvar, uvar)
     Xpress.lpoptimize(prob)
     # x, slack, dual, dj = Float64[], Float64[], Float64[], Float64[]
     # Xpress.Lib.XPRSgetsol(prob, x, slack, dual, dj)
     # if QM.meta.nnzh > 0
end

end

end
