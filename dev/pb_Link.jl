using DataStructures

struct Link
  orig::String
  dest::String
end

import Base.isless

function isless(l1::Link, l2::Link)
    return isless((l1.orig, l1.dest), (l2.orig, l2.dest))
end

d = SortedDict{Link, SortedDict{String, Any}}()

l_test = Link("o", "d")
d[l_test] = SortedDict{String, Any}()


d[l_test]["elem1"] = 1

haskey(d, Link("o", "d"))

d[Link("o", "d")]["elem2"] = "test"
