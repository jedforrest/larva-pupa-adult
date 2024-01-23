
using HomotopyContinuation
@var x y
f = System([x^2 + y^2 - 4, x - y])
result = solve(f)
open("julia_output.txt", "w") do io
    println(io, real_solutions(result))
end
