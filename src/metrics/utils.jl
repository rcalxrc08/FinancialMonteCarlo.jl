# A multivariate model knows neither its name nor those of its components, so the positions priced on it are given by
# subset of its components, in a fixed order: the nonempty subsets of 1:n by size, then lexicographically. For n = 3:
# [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]. A subset is the set bits of a number in 1:2^n - 1.
compute_indices(n::Integer) = sort!([findall(i -> isodd(mask >> (i - 1)), 1:n) for mask = 1:(2^n-1)]; by = subset -> (length(subset), subset))
