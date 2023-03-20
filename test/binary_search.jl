struct ValidityCheck
    t::Int
end

(v::ValidityCheck)(t) = t ≤ v.t

for i ∈ 0:10
    v,t = binary_search_max(identity, ValidityCheck(i), 10)
    @test v == t == i
end
