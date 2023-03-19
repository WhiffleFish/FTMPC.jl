struct ValidityCheck
    t::Int
end

(v::ValidityCheck)(t) = t ≤ v.t

binary_search_max(identity, ValidityCheck(5), 20)
