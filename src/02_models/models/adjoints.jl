# Adjoints (reverse-mode derivatives) of the kernels of the paths, in closed form and on whole arrays. They are shared by
# the reverse rules of the Mooncake and Enzyme extensions; the rules only fetch the arrays and call these functions.
#
# Why: without them both backends differentiate the broadcasts of the kernels element by element. Mooncake pushes the
# data of every element on its stack, and Enzyme caches every value a kernel overwrites (the normals are redrawn into the
# same buffer at every step). Here each step is one pass over the paths.
#
# Notation and conventions, for a kernel y = f(x, p) on the paths (x arrays, p scalars) and a scalar loss L:
# - ȳ = ∂L/∂y is the adjoint of the output, given: it is the shadow (tangent buffer) of the array the kernel wrote,
#   called dy (or dx, dw, ... for the array of that name).
# - The adjoint of an array input x accumulates x̄ += (∂y/∂x)ᵀ ȳ, element by element: with one output per path,
#   x̄_i += ȳ_i ∂y_i/∂x_i. An input that is constant for the backend (Enzyme's Const) has no shadow and is passed as
#   `nothing`: accumulate! then does nothing.
# - The adjoint of a scalar parameter p is the sum over the paths p̄ = Σ_i ȳ_i ∂y_i/∂p: it is returned, and the rule hands
#   it back to the backend.
# - When the kernel writes over an array, the old values reach the loss only through the kernel. So the shadow of that
#   array must end up as the adjoint of the old values: zero if the old values were not read (a new column of the paths),
#   or the propagated adjoint if they were (the variance of Heston, the log paths exponentiated in place).
# - The values the formulas need (the normals of the step, the variance before the step, ...) are passed in by the rules,
#   which saved them in the forward pass, since the arrays may have been overwritten since then.

# x̄[i] += v, for the shadow of an array input; nothing to do for a constant input (no shadow)
@inline accumulate!(s, i, v) = (s[i] += v; nothing)
@inline accumulate!(::Nothing, i, v) = nothing

# bm_step!: one step of a Brownian motion on the paths, y_i = x_i + drift + s z_i (y the new column, x the previous one,
# z the normals, s = ±σ√dt the signed standard deviation of the step).
# Partial derivatives: ∂y_i/∂x_i = 1, ∂y_i/∂drift = 1, ∂y_i/∂s = z_i, ∂y_i/∂z_i = s. Hence
#   x̄_i += ȳ_i,   z̄_i += s ȳ_i,   drift̄ = Σ ȳ_i,   s̄ = Σ ȳ_i z_i = ȳ ⋅ z.
# The new column overwrote values that were never read, so its adjoint is reset to zero.
function bm_step_adjoint!(dx, dz, dy, z, s)
    ds = dot(dy, z)
    ddrift = sum(dy)
    @inbounds for i in eachindex(dy, z)
        accumulate!(dx, i, dy[i])
        accumulate!(dz, i, s * dy[i])
        dy[i] = zero(eltype(dy))
    end
    return ddrift, ds
end

# exp_paths!: the prices from the log paths, in place, y = S0 exp(x).
# Partial derivatives: ∂y/∂x = S0 exp(x) = y, ∂y/∂S0 = exp(x) = y / S0. Hence
#   x̄ = ȳ y   (the shadow of X turns from the adjoint of the prices into the adjoint of the log paths),
#   S̄0 = Σ ȳ exp(x) = (ȳ ⋅ y) / S0, from the old values x_old when S0 = 0 (where y / S0 is undefined).
function exp_paths_adjoint!(dx, x_new, x_old, s0)
    ds0 = iszero(s0) ? sum(i -> dx[i] * exp(x_old[i]), eachindex(dx, x_old)) : dot(dx, x_new) / s0
    dx .*= x_new
    return ds0
end

# affine_carry!: the paths of Bachelier, in place, y_ij = (a + x_ij) c_j, with x the arithmetic Brownian motion, a = S0 and
# c_j = e^{∫₀^{t_j} (r - d)} the carry factor of the column j.
# Partial derivatives: ∂y_ij/∂x_ij = c_j, ∂y_ij/∂a = c_j, ∂y_ij/∂c_j = a + x_ij. Hence
#   x̄_ij = ȳ_ij c_j   (the shadow of X turns into the adjoint of the old values),
#   ā = Σ_ij ȳ_ij c_j,   c̄_j += Σ_i ȳ_ij (a + x_ij) with the old values x_old.
@muladd function affine_carry_adjoint!(dx, dc, x_old, a, c)
    da = zero(eltype(dx))
    @inbounds for j in axes(dx, 2)
        dcj = zero(eltype(dx))                 # c̄_j, summed over the paths
        for i in axes(dx, 1)
            x̄ = dx[i, j]                       # ȳ_ij
            dcj = dcj + x̄ * (a + x_old[i, j])
            da = da + x̄ * c[j]
            dx[i, j] = x̄ * c[j]                # x̄_ij
        end
        accumulate!(dc, j, dcj)
    end
    return da
end

# add_carry!: the deterministic shift of the shifted models, in place, y_ij = x_ij + b c_j, with b = α S0 and c_j the carry
# factor of the column j.
# Partial derivatives: ∂y_ij/∂x_ij = 1, ∂y_ij/∂b = c_j, ∂y_ij/∂c_j = b. Hence
#   x̄_ij = ȳ_ij   (the shadow of X is unchanged),   b̄ = Σ_j c_j Σ_i ȳ_ij,   c̄_j += b Σ_i ȳ_ij.
@muladd function add_carry_adjoint!(dx, dc, b, c)
    db = zero(eltype(dx))
    @inbounds for j in axes(dx, 2)
        column = zero(eltype(dx))              # Σ_i ȳ_ij
        for i in axes(dx, 1)
            column += dx[i, j]
        end
        db = db + column * c[j]
        accumulate!(dc, j, b * column)
    end
    return db
end

# heston_step!: one Euler step of Heston on the paths of one side. For each path, with v the variance before the step,
# n1, n2 the two normals, a = ±√dt and fl the floor of the variance:
#   y   = x + drift - v dt / 2 + √v a n1                          (the log price, the new column)
#   e2  = ρ n1 + ρo n2                                             (the normal of the variance, correlated)
#   raw = v + κ (θ - v) dt + σ √v a e2,   v_new = max(raw, fl)     (the variance, in place)
# The adjoints of the outputs are ȳ (dy) and v̄_new (dw). The floor passes v̄_new to raw where raw > fl, and to fl
# elsewhere: g = v̄_new where the variance stayed above the floor, 0 where the floor bound.
# Partial derivatives of y:   ∂/∂x = 1, ∂/∂drift = 1, ∂/∂v = -dt / 2 + a n1 / (2√v), ∂/∂dt = -v / 2,
#                             ∂/∂a = √v n1, ∂/∂n1 = √v a.
# Partial derivatives of raw: ∂/∂v = 1 - κ dt + σ a e2 / (2√v), ∂/∂dt = κ (θ - v), ∂/∂a = σ √v e2, ∂/∂κ = (θ - v) dt,
#                             ∂/∂θ = κ dt, ∂/∂σ = √v a e2, ∂/∂e2 = σ √v a, with ∂e2/∂n1 = ρ, ∂e2/∂n2 = ρo,
#                             ∂e2/∂ρ = n1, ∂e2/∂ρo = n2.
# Each adjoint is ȳ times the partial of y plus g times the partial of raw. The shadow of v (dw) turns from the adjoint of
# the new variance into the adjoint of the old one, w_old the variance before the step.
@muladd function heston_step_adjoint!(dx, dw, dz1, dz2, dy, w_old, z1, z2, drift, dt, a, κ, θ, σ, ρ, ρo, fl)
    P = eltype(dy)
    ddrift = ddt = da = dκ = dθ = dσ = dρ = dρo = dfl = zero(P)
    @inbounds for i in eachindex(dy, dw, w_old, z1, z2)
        ȳ, w̄ = dy[i], dw[i]                       # adjoints of the new log price and of the new variance
        vi, n1, n2 = w_old[i], z1[i], z2[i]       # the variance before the step and the two normals
        sv = sqrt(vi)
        e2 = n1 * ρ + n2 * ρo
        # same test as the kernel (the same scalar function, so the same rounding): did the floor bind?
        alive = fl < heston_variance_step(vi, n1, n2, dt, a, κ, θ, σ, ρ, ρo)
        g = alive ? w̄ : zero(P)                   # adjoint of raw
        dfl += alive ? zero(P) : w̄                # adjoint of the floor, where it bound
        accumulate!(dx, i, ȳ)
        # adjoint of the old variance: through the log price and through raw
        dw[i] = ȳ * (a * n1 / (2 * sv) - dt / 2) + g * (1 - κ * dt + σ * a * e2 / (2 * sv))
        # the normals: n1 enters y directly and e2 with ρ, n2 only e2 with ρo
        accumulate!(dz1, i, sv * a * (ȳ + g * σ * ρ))
        accumulate!(dz2, i, g * σ * sv * a * ρo)
        dy[i] = zero(P)                           # the new column overwrote values never read
        # the scalars, summed over the paths
        ddrift += ȳ
        ddt = ddt - ȳ * vi / 2 + g * κ * (θ - vi)
        da = da + sv * (ȳ * n1 + g * σ * e2)
        dκ = dκ + g * (θ - vi) * dt
        dθ = dθ + g * κ * dt
        dσ = dσ + g * sv * a * e2
        dρ = dρ + g * σ * sv * a * n1
        dρo = dρo + g * σ * sv * a * n2
    end
    return ddrift, ddt, da, dκ, dθ, dσ, dρ, dρo, dfl
end

# lnm_step!: one Euler step of the lognormal mixture on the paths of one side. For each path, with x the log price before
# the step and z = sign Z its normal:
#   y_new = x + drift - σ² h / 2 + √(σ² h) z,   σ² = σ²(t, x - carry) the local variance (the constant σ2_0 at t = 0).
# The local variance is a weighted mean of the η_k², with weights w_k ∝ λ_k p_k (p_k the density of the component k).
# At the log-forward-moneyness ψ = x - carry, written with u_k = ψ + η_k² t / 2 and logw_k = log(λ_k / η_k):
#   e_k = logw_k - u_k² / (2 η_k² t),   w_k = exp(e_k - max e),   W = Σ w_k,   σ² = Σ w_k η_k² / W.
# lnm_variance_partials returns σ² and its partial derivatives:
#   ∂σ²/∂e_k = w_k (η_k² - σ²) / W                     (the derivative of a weighted mean in the log weights)
#   ∂e_k/∂ψ = -u_k / (η_k² t),  ∂e_k/∂t = u_k (u_k / (η_k² t) - 1) / (2t),  ∂e_k/∂η_k = u_k (u_k / (η_k² t) - 1) / η_k,
#   ∂e_k/∂logw_k = 1, and σ² depends on η_k also directly, through η_k²: + 2 η_k w_k / W.
# The chain rule gives ∂σ²/∂ψ (dy), ∂σ²/∂t (dt), ∂σ²/∂η_k (dη) and ∂σ²/∂logw_k (de).
@muladd function lnm_variance_partials(y, t, η::NTuple{N}, logw::NTuple{N}) where {N}
    u = map(η_k -> y + η_k^2 * (t / 2), η)
    e = map((u_k, η_k, l_k) -> l_k - u_k^2 / (2 * η_k^2 * t), u, η, logw)
    e_max = maximum(e)
    w = map(e_k -> exp(e_k - e_max), e)
    W = sum(w)
    σ2 = sum(map((w_k, η_k) -> w_k * η_k^2, w, η)) / W
    de = map((w_k, η_k) -> w_k * (η_k^2 - σ2) / W, w, η)                                   # ∂σ²/∂e_k
    dy = -sum(map((d, u_k, η_k) -> d * u_k / (η_k^2 * t), de, u, η))                         # ∂σ²/∂ψ
    dt = sum(map((d, u_k, η_k) -> d * u_k * (u_k / (η_k^2 * t) - 1) / (2 * t), de, u, η))     # ∂σ²/∂t
    dη = map((d, w_k, u_k, η_k) -> 2 * η_k * w_k / W + d * u_k * (u_k / (η_k^2 * t) - 1) / η_k, de, w, u, η) # ∂σ²/∂η_k
    return σ2, dy, dt, dη, de
end
# The local variance and its partials at one path; at t = 0 the constant σ2_0, which depends on nothing on the path
function lnm_variance_and_partials(x, t, carry, σ2_0::P, η::NTuple{N, P}, logw::NTuple{N, P}) where {N, P}
    iszero(t) && return σ2_0, zero(P), zero(P), ntuple(_ -> zero(P), Val(N)), ntuple(_ -> zero(P), Val(N))
    return lnm_variance_partials(x - carry, t, η, logw)
end
# Partial derivatives of y_new:
#   ∂/∂σ² = -h / 2 + z √(σ² h) / (2σ²)    (q / ȳ below),   ∂/∂h = -σ² / 2 + z √(σ² h) / (2h),
#   ∂/∂x = 1 + ∂/∂σ² ∂σ²/∂ψ,   ∂/∂carry = -∂/∂σ² ∂σ²/∂ψ,   ∂/∂drift = 1,   ∂/∂Z = sign √(σ² h),
#   and through σ²: ∂/∂t, ∂/∂η_k, ∂/∂logw_k (t > 0) or ∂/∂σ2_0 (t = 0).
@muladd function lnm_step_adjoint!(dx, dz, dy, x, z, sign, t, h, drift, carry, σ2_0::P, η::NTuple{N, P}, logw::NTuple{N, P}) where {N, P}
    dt = dh = ddrift = dcarry = dσ2_0 = zero(P)
    dη = dlogw = ntuple(_ -> zero(P), Val(N))
    @inbounds for i in eachindex(dy, x, z)
        ȳ = dy[i]
        σ2, py, pt, pη, pl = lnm_variance_and_partials(x[i], t, carry, σ2_0, η, logw)
        zi = sign * z[i]
        sq = sqrt(σ2 * h)
        q = ȳ * (zi * sq / (2 * σ2) - h / 2)       # ȳ ∂y/∂σ², the adjoint of the local variance
        accumulate!(dx, i, ȳ + q * py)
        accumulate!(dz, i, ȳ * sign * sq)
        ddrift += ȳ
        dh = dh + ȳ * (zi * sq / (2 * h) - σ2 / 2)
        if iszero(t)
            dσ2_0 += q
        else
            dcarry = dcarry - q * py
            dt = dt + q * pt
            dη = map((acc, p) -> acc + q * p, dη, pη)
            dlogw = map((acc, p) -> acc + q * p, dlogw, pl)
        end
        dy[i] = zero(P)                            # the new column overwrote values never read
    end
    return dt, dh, ddrift, dcarry, dσ2_0, dη, dlogw
end

# add_jump_sizes!: the jumps of Kou and Merton, X[path[k], j] += sizes[k] for every grid point j ≥ step[k] (a jump moves
# the rest of its path). Partial derivatives: ∂X_ij/∂X_ij = 1, ∂X_ij/∂sizes[k] = 1 on the row path[k] from step[k] on.
# Hence the shadow of X is unchanged and sizes̄[k] += Σ_{j ≥ step[k]} X̄[path[k], j]. The grid points and the paths carry
# no derivative (integers).
function add_jump_sizes_adjoint!(dsizes, dX, path, step)
    @inbounds for k in eachindex(path, step, dsizes)
        i = path[k]
        acc = zero(eltype(dsizes))
        for j = step[k]:size(dX, 2)
            acc += dX[i, j]
        end
        dsizes[k] += acc
    end
    return nothing
end
