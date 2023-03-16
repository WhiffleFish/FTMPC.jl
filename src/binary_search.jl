function binary_search_max(f, isvalid, T)
    L = 0
    R = T-1
    val = nothing
    while L ≤ R
        t = (L + R) ÷ 2
        val = f(t)
        if isvalid(val)
            L = t + 1
        else
            R = t - 1
        end
    end
    return val, L
end
