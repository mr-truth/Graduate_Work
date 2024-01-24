using QuadGK

na_array = [(10^i) for i in -4:0.2:-2]
pump = 0.0625
t = 0
compensation = 0.025
E_array = [i for i in -2:0.001:6]
PL = zeros(length(E_array))
E_max = zeros(length(na_array))
x_max = zeros(length(na_array))

for k=1:length(na_array)
na = na_array[k]
E_md = 2.515 * 2 * compensation * (na^(2/3))
println(na)
	for i=1:length(E_array)
		E = E_array[i]
		
		function integrate()
			func(energy_1) = (energy_1^-5) * (E_md/((4*(E_md/energy_1)^2) + (E - energy_1)^2)) * exp(-4 / energy_1 - 32 * π / 3 * na / energy_1 ^ 3 - t*exp(-4/energy_1)) 
			result = quadgk(func, 0, 1000, rtol=1e-12)[1]
			end
		
		function integrate_q(t)
			func(energy) = (exp(-t * exp(-4/energy)) - 1)/energy^4
			result = quadgk(func, 0, 1000, rtol=1e-12)[1]
			end
		
		a = integrate() * exp(integrate_q(t) * 4 * π *na)
		PL[i] = (64 * na * pump)/(32 * π * na) * a
		
		println(PL[i])
	end
max = maximum(PL)
E_max[k] = max
for (x, y) in zip(E_array, PL)
	if y == max
		E_max[k] = y
		x_max[k] = x
	end
end
println("--------------")
end

for i = 1:length(E_max)
	e_max = E_max[i]
	println(e_max)
end

println("************")

for i = 1:length(x_max)
	x = x_max[i]
	println(x)
end

println("************")

for i = 1: length(na_array)
	n = na_array[i]
	println(n)
end