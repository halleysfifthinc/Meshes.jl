# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type AbstractConnectivity end

"""
    Connectivity{PL,N}

A connectivity list of `N` indices representing a [`Polytope`](@ref)
of type `PL`. Indices are taken from a global vector of [`Point`](@ref).

Connectivity objects are constructed with the [`connect`](@ref) function.
"""
struct Connectivity{PL<:Polytope,N} <: AbstractConnectivity
  indices::NTuple{N,Int}

  function Connectivity{PL,N}(indices) where {PL,N}
    assertion(nvertices(PL) == N, lazy"cannot create a $PL with $N vertices")
    if PL <: Ngon
      assertion(N ≥ 3, "Ngon requires 3 or more vertices")
    end
    new(indices)
  end
end

struct PolyConnectivity <: AbstractConnectivity
    indices::Vector{Int}
end

"""
    pltype(connectivity)

Return the polytope type that the connectivity represents.
"""
pltype(c::AbstractConnectivity) = pltype(typeof(c))
pltype(::Type{Connectivity{PL,N}}) where {PL,N} = PL
pltype(::Type{PolyConnectivity}) = PolyArea

"""
    paramdim(connectivity)

Return the parametric dimension of the `connectivity`.
"""
paramdim(c::AbstractConnectivity) = paramdim(typeof(c))
paramdim(::Type{Connectivity{PL,N}}) where {PL,N} = paramdim(PL)
paramdim(::Type{PolyConnectivity}) = paramdim(PolyArea)

"""
    indices(connectivity)

Return the list of indices of the `connectivity`.
"""
indices(c::AbstractConnectivity) = c.indices

"""
    connect(indices, [PL])

Connect a list of `indices` from a global vector of [`Point`](@ref)
into a [`Polytope`](@ref) of type `PL`.

The type `PL` can be a [`Ngon`](@ref) in which case the length of
the indices is used to identify the actual polytope type.

Finally, the type `PL` can be ommitted. In this case, the indices are
assumed to be connected as a [`Ngon`](@ref) or as a [`Segment`](@ref).

## Examples

Connect indices into a Triangle:

```julia
connect((1,2,3), Triangle)
```

Connect indices into N-gons, a `Triangle` and a `Quadrangle`:

```julia
connect.([(1,2,3), (2,3,4,5)], Ngon)
```

Connect indices into N-gon or segment:

```julia
connect((1,2)) # Segment
connect((1,2,3)) # Triangle
connect((1,2,3,4)) # Quadrangle
```
"""
connect(indices::Tuple, PL::Type{<:Polytope}) = Connectivity{PL,length(indices)}(indices)

connect(indices, ::Type{PolyArea}) = PolyConnectivity(indices)
connect(indices::Tuple, ::Type{PolyArea}) = PolyConnectivity(collect(indices))

function connect(indices::Tuple, ::Type{Ngon})
  N = length(indices)
  Connectivity{Ngon{N},N}(indices)
end

function connect(indices::Tuple)
  N = length(indices)
  N > 2 ? connect(indices, Ngon) : connect(indices, Segment)
end

function connect(indices::AbstractVector{<:Integer})::Connectivity
  connect(ntuple(i -> Int(indices[i]), length(indices)), Ngon)
end

"""
    materialize(connec, points)

Materialize a face using the `connec` list and a global vector of `points`.
"""
materialize(connec::Connectivity{PL,N}, points::AbstractVector{P}) where {PL<:Polytope,N,P<:Point} =
  PL(ntuple(i -> @inbounds(points[connec.indices[i]]), N))

function materialize(connec::PolyConnectivity, points::AbstractVector{P}) where {P<:Point}
  PolyArea(map(i -> points[i], indices(connec)))
end


function Base.show(io::IO, c::Connectivity{PL}) where {PL}
  name = prettyname(PL)
  inds = c.indices
  print(io, "$name$inds")
end
