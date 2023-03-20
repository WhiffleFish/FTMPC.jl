function binary_search_max(f, isvalid, T)
    best_t = 0
    _min = 0
    _max = T
    best_val = nothing
    while _min < _max
        t = (_min + _max) ÷ 2
        val = f(t)
        if isvalid(val) # consensus horizon too low
            _min = t + 1
            best_t = t
            best_val = val
        else            # consensus horizon too high
            _max = t - 1
        end
    end
    max_val = f(_max)
    return if isvalid(max_val)
        max_val, _max
    else
        best_val, best_t
    end
end
