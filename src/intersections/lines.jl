# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of three types:
#
# 1. intersect at one point
# 2. overlap at more than one point
# 3. do not overlap nor intersect
function intersection(f, line₁::Line, line₂::Line)
  a, b = line₁(0), line₁(1)
  c, d = line₂(0), line₂(1)

  λ₁, _, r, rₐ = intersectparameters(a, b, c, d)

  if r == rₐ == 2
    return @IT Crossing (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT Overlapping line₁ f
  else
    return @IT NotIntersecting nothing f
  end
end

function intersection(f, line::Line, poly::PolyArea)
  segs = Iterators.flatten(segments.(rings(poly)))
  ints = Intersection[]

  for seg in segs
    int = intersection(identity, seg, line)
    if type(int) !== NotIntersecting
      push!(ints, int)
    end
  end
  isempty(ints) && return @IT NotIntersecting nothing f

  # delete duplicate intersections from end/begining of segments on either side of an
  # overlapping segment. depending on ring orientation, the duplicate points in this
  # specific case can be hidden from the simple de-dupping by the upcoming sorting.
  # duplicate end/beginning without an overlapping segment in-between will still be adjacent
  # after sorting, so don't need to be handled here
  i = 1
  while i < lastindex(ints)
    int = ints[i]
    if type(int) === Touching
      if i+1 < lastindex(ints)
        int1 = ints[i+1]
        if type(int1) === Overlapping
          if i+2 ≤ lastindex(ints)
            int2 = ints[i+2]
            if type(int2) === Touching
              pt = get(int)
              pt1_1, pt1_2 = (vertices(get(int1))...,)
              pt2 = get(int2)

              if pt ≈ pt1_1 && pt1_2 ≈ pt2
                  deleteat!(ints, (i, i+2))
                  continue
              end
            end
          end
        end
      end
    end
    i += 1
  end

  # Flip segments as needed so the first vertex is closer to ray origin
  uv = unormalize(line(1) - line(0))
  for (i,int) in enumerate(ints)
    if type(int) === Overlapping
      v1, v2 = (vertices(get(int))...,)
      if udot(to(v1), uv) > udot(to(v2), uv)
        ints[i] = @IT Overlapping Segment(v2, v1) identity
      end
    end
  end
  sort!(ints; by=function(int)
        pt = if type(int) === Overlapping
          vertices(get(int))[1]
        else
          get(int)
        end
        udot(to(pt), uv)
    end)


  PTT = eltype(vertices(poly))
  RT = Rope{manifold(poly),crs(poly),Vector{PTT}}
  set = [RT(type(first(ints)) === Overlapping ? vertices(get(first(ints))) : [get(first(ints))])]

  for int in @view(ints[2:end])
    pt_prev = last(last(set).vertices)
    pt = type(int) === Overlapping ? first(vertices(get(int))) : get(int)
    seg = Segment(pt_prev, pt)

    if pt_prev ≈ pt || seg(0.5) ∈ poly
      pt_prev ≈ pt || push!(last(set).vertices, pt)
      if type(int) === Overlapping
        push!(last(set).vertices, last(vertices(get(int))))
      end
    else
      if type(int) === Overlapping
        p2 = last(vertices(get(int)))
        push!(set, RT([pt, p2]))
      else
        push!(set, RT([pt]))
      end
    end
  end

  return @IT Intersecting GeometrySet(set) f
end

