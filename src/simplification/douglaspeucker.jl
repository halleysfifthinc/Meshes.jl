# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DouglasPeuckerSimplification(τ)

Douglas-Peucker's simplification algorithm with tolerance `τ` in length units
(default to meter).

The higher is the tolerance, the more aggressive is the simplification.

## References

* Douglas, D. and Peucker, T. 1973. [Algorithms for the Reduction of
  the Number of Points Required to Represent a Digitized Line or its
  Caricature](https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
struct DouglasPeuckerSimplification{ℒ<:Len} <: SimplificationMethod
  τ::ℒ
  DouglasPeuckerSimplification(τ::ℒ) where {ℒ<:Len} = new{float(ℒ)}(τ)
end

DouglasPeuckerSimplification(τ) = DouglasPeuckerSimplification(addunit(τ, u"m"))

function simplify(chain::Chain, method::DouglasPeuckerSimplification)
  verts = _douglaspeucker(vertices(chain), method.τ)
  isclosed(chain) ? Ring(verts) : Rope(verts)
end

# simplify chain assuming it is open
function _douglaspeucker(v::AbstractVector{P}, τ) where {P<:Point}
  n = length(v)

  if n ≤ 2
      return copy(v)
  end

  keep = BitSet()
  push!(keep, 1)
  push!(keep, n)

  stack = Vector{Tuple{Int, Int}}()
  push!(stack, (1, n))

  while !isempty(stack)
    start_idx, end_idx = pop!(stack)

    # find vertex with maximum distance to reference line
    l = Line(v[start_idx], v[end_idx])
    imax, dmax = 0, zero(lentype(P))
    for i in (start_idx + 1):(end_idx - 1)
      d = evaluate(Euclidean(), v[i], l)
      if d > dmax
        imax = i
        dmax = d
      end
    end

    # if maximum distance exceeds tolerance, split segment
    if dmax ≥ τ
      push!(keep, imax)
      # only push sub-chains if they are splittable (i.e. have a point between endpoints)
      imax > start_idx+1 && push!(stack, (start_idx, imax))
      imax < end_idx-1 && push!(stack, (imax, end_idx))
    end
  end

  # BitSet's are already sorted, so the indices will collect elements of `v` in the original
  # order
  return v[collect(keep)]
end

