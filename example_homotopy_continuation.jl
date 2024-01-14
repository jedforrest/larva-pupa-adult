using HomotopyContinuation
@var x y
f = System([
		   x^2 + y^2 - 4,
		   x - y,
])
result = HomotopyContinuation.solve(f)
sols = real_solutions(result)
if length(sols) > 0
	open("hc_output.txt", "w") do io
		println(io, real_solutions(result))
	end
end
