function binary_search_max(f, isvalid, T)
    L = 0
    R = T-1
    best_val = nothing
    while L < R
        t = (L + R) ÷ 2
        val = f(t)
        if isvalid(val) # consensus horizon too high
            L = t + 1
            best_val = val
        else            # consensus horizon too low
            R = t - 1
        end
    end
    return best_val, L
end
