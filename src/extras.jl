function basis(n,i)
    e_i = zeros(n)
    e_i[i] = 1
    return e_i
end

function blkdiag(mat::AbstractMatrix, n::Int)
    return cat(Iterators.repeated(mat, n)..., dims=(1,2))
end

blkdiag(mats) = cat(mats..., dims=(1,2))

convert_kwargs(kwargs::Base.Pairs) = Tuple(string(a)=>b for (a,b) ∈ kwargs)
