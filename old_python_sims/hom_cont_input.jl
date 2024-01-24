
using HomotopyContinuation
@var b c_ea c_el c_pa mu_a
fL = [-157.305487948269*b*c_ea - 393.263719870674*b*c_el + 8.10461575772587*b - 10.1192217506699, -2868.81807229302*b*c_ea - 283.603042721896*b*c_el + 64.330653583117*b - 183.214836377097, -424.668177446337*b*c_ea - 346.801545418276*b*c_el + 10.9315442189046*b - 12.1054055320091]
fA = [-24.1341011269039*c_pa - 100*mu_a + 2.63485429691197; -9725.73748664558*c_pa - 102.362111341413*mu_a - 4.45779910242385; -61.0753833435726*c_pa - 224.351683760664*mu_a + 6.9379620128420]
f = System(fA)
result = HomotopyContinuation.solve(f)
real_solutions(result)
open("julia_output.txt", "w") do io
    println(io, real_solutions(result))
end
