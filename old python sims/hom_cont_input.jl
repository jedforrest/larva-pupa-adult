
using HomotopyContinuation
@var b c_ea c_el
f = System([-157.305487948269*b*c_ea - 393.263719870674*b*c_el + 8.10461575772587*b - 10.1192217506699, -2868.81807229302*b*c_ea - 283.603042721896*b*c_el + 64.330653583117*b - 183.214836377097, -424.668177446337*b*c_ea - 346.801545418276*b*c_el + 10.9315442189046*b - 12.1054055320091])
result = solve(f)
open("julia_output.txt", "w") do io
    println(io, real_solutions(result))
end
