################################################################################################################
#
# bounded_diff.jl = a julia-language script for modeling transient 1-D gaseous diffusion through
# partially water-saturated porous media, taking into account soil water-solid partitioning onto organic carbon
#
################################################################################################################

type Solute
	D::Float64 				# diffusion coefficient
	H::Float64 				# dimensionless Henry's law coefficient
	Koc::Float64 			# organic-carbon partitioning coefficient
end


type Bar
	L::Float64 				# length of 1-D domain
	area::Float64 			# cross-sectional area
	phi::Float64 			# porosity
	rho_b::Float64 			# soil bulk density
	theta::Float64 			# fractional pore water content
	foc::Float64 			# organic carbon fraction
	C1::Float64 			# boundary concentration; C(0, t) = C1; for other boundary, C(L, t) = 0.0
end


type Knobs
	gamma::Float64 			# time-stepping weighting factor for implicit solution scheme
	dt_init::Float64 		# initial, minimum, and maximum tim estep size
	dt_min::Float64
	dt_max::Float64
	du_max::Float64 		# maximum change in concentration, per time step
	dt_decrease::Float64	# time step reduction and increase factors
	dt_increase::Float64
end


type Node
	x::Float64 					# node point location
	V::Float64 					# cell volume
	phi::Float64 				# porosity
	theta::Float64 				# water content
	sigma::Float64 				# cell total conductance (summed over all connections)
	S::Float64 					# source term (mass flux)
	u::Float64 					# concentration
end


type Connection
	node_1::Int64 				# index numbers for connecting nodes
	node_2::Int64
	dx::Float64					# inter-node distance (for computing gradient)
	area::Float64 				# cell connection interfacial area
	diff_conduct::Float64 		# connection conductance incorporates dx and area; these will be used in subsequent expansions of the model
end


function GetKnobs()
	# read numerical model "knobs" from file
	data = readdlm("knobs.txt", '\t', header=false)
	gamma = Float64(data[1, 2])
	dt_init = Float64(data[2, 2])
	dt_min = Float64(data[3, 2])
	dt_max = Float64(data[4, 2])
	du_max = Float64(data[5, 2])
	dt_decrease = Float64(data[6, 2])
	dt_increase = Float64(data[7, 2])
	knobs = Knobs(gamma, dt_init, dt_min, dt_max, du_max, dt_decrease, dt_increase)
	return knobs
end


function GetParams()
	# read basic model parameters from file
	data = readdlm("parameters.txt", '\t', header=false)
	D = Float64(data[1, 2])
	H = Float64(data[2, 2])
	Koc = Float64(data[3, 2])
	L = Float64(data[4, 2])
	area = Float64(data[5, 2])
	C1 = Float64(data[6, 2])
	phi = Float64(data[7, 2])
	sat = Float64(data[8, 2])
	rho_b = Float64(data[9, 2])
	foc = Float64(data[10, 2])
	# assign parameter values to appropriate types
	theta = sat * phi
	solute = Solute(D, H, Koc)
	bar = Bar(L, area, phi, rho_b, theta, foc, C1)
	return solute, bar
end


function R_coeff(solute::Solute, bar::Bar)
	# effective species retardation coefficient with respect to gas phase
	return 1.0 + bar.theta/(bar.phi-bar.theta) * (1/solute.H) + bar.rho_b/(bar.phi-bar.theta) * solute.Koc * bar.foc * (1/solute.H)
end


function C_a(x::Float64, t::Float64, D::Float64, R::Float64, bar::Bar)
	# analytical solution for C(x, t) for a semi-infinite half-space (good approximation for 1-D domain at early time)
	u = bar.C1 * erfc(x/(2*sqrt(D/R*t)))
	return u 		# returns a scalar (at x)
end


function LHS_matrix(connection::Array{Connection,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64)
	# fill out the LHS of the equation matrix (off-diagonal terms, listed last, will not change in a static materials problem)
	row_index = Int64[] 					# indexing system for sparse matrix
	col_index = Int64[]
	data = Float64[]
	for (i, nd) in enumerate(node)         	# diagonal terms
		push!(row_index, i)
		push!(col_index, i)
		push!(data, nd.V/dt + knobs.gamma*nd.sigma)
	end
	for connect in connection 				# off-diagonal terms
		push!(row_index, connect.node_1)
		push!(col_index, connect.node_2)
		push!(data, -knobs.gamma * connect.diff_conduct)
		push!(row_index, connect.node_2)
		push!(col_index, connect.node_1)
		push!(data, -knobs.gamma * connect.diff_conduct)
	end
	return data, row_index, col_index
end


function Diagonal(data::Array{Float64,1}, node::Array{Node,1}, knobs::Knobs, dt::Float64)
	# update equation matrix diagonal
	for (i, nd) in enumerate(node)
		data[i] = nd.V/dt + knobs.gamma*nd.sigma
	end
	return data
end


function RHS_vector(connection::Array{Connection,1}, node::Array{Node,1})
	# construct explicit matrix (run for each time step)
	b = Float64[]
	for (i, nd) in enumerate(node)
		push!(b, nd.S - nd.sigma * nd.u)
	end	
	for connect in connection
        b[connect.node_1] += connect.diff_conduct * node[connect.node_2].u
		b[connect.node_2] += connect.diff_conduct * node[connect.node_1].u
	end	
	return b
end


function C_n(x::Array{Float64,1}, t_end::Float64, D_scalar::Float64, R_scalar::Float64, bar::Bar, knobs::Knobs, num_pts::Int64)

	# numerical solution for C(x, t) inside bar; simplified to assume homogeneity

	# initialize arrays
	dx = bar.L/num_pts
	num_nodes = num_pts + 2 										# interior nodes, plus boundaries	
	V = zeros(Float64, num_nodes) + dx * bar.area 					# cell volume
	phi = zeros(Float64, num_nodes) + bar.phi
	theta = zeros(Float64, num_nodes) + bar.theta
	S = zeros(Float64, num_nodes) 									# source term (not used for this simple application)
	u = zeros(Float64, num_nodes) 									# initial conditions
	
	# designate boundary nodes by assigning large volumes
	V[1] = 1e+20
	u[1] = bar.C1
	V[num_nodes] = 1e+20
	u[num_nodes] = 0.0
	
	# set up connections (this feature provides flexibility for irregular column delineation)
	connection = Connection[]
	push!(connection, Connection(1, 2, dx, bar.area, (bar.phi-bar.theta)*D_scalar*bar.area/(R_scalar*dx)))					# first boundary
	push!(connection, Connection(num_nodes-1, num_nodes, dx, bar.area, (bar.phi-bar.theta)*D_scalar*bar.area/(R_scalar*dx))) 	# second boundary
	for i = 2:num_nodes-1 		# interior nodes
		push!(connection, Connection(i-1, i, dx, bar.area, (bar.phi-bar.theta)*D_scalar*bar.area/(R_scalar*dx)))
	end
	
	# summarize conductances, per node
	sigma = zeros(Float64, num_nodes)
	for connect in connection
		sigma[connect.node_1] += connect.diff_conduct
		sigma[connect.node_2] += connect.diff_conduct		
	end

	# set up node data type array
	node = Node[]
	push!(node, Node(x[1]-dx, V[1], phi[1], theta[1], sigma[1], S[1], u[1]))
	for i = 2:num_nodes-1
		push!(node, Node(x[i-1], V[i], phi[i], theta[i], sigma[i], S[i], u[i]))
	end
	push!(node, Node(x[num_pts]+dx, V[num_pts], phi[num_pts], theta[num_pts], sigma[num_pts], S[num_pts], u[num_pts]))	
	
	t = 0.
	dt = knobs.dt_init
	
	# set up left-hand-side of flux balance equations matrix
	data, row_index, col_index = LHS_matrix(connection, node, knobs, dt)
	
    while (t < t_end)

        complete = false

        while complete == false

            data = Diagonal(data, node, knobs, dt) 							# modify diagonal according to time step size
 			A = sparse(row_index, col_index, data, num_nodes, num_nodes) 	# update sparse equation matrix
            b = RHS_vector(connection, node) 								# construct explicit vector
			global du = \(A, b)	 											# solve equation set	
			
            # check maximum concentration change criterion at this time step size
			sum_complete = 0
			for i = 1:num_nodes
				sum_complete += (abs(du[i]) > knobs.du_max)
			end
			complete = 1 - sign(sum_complete)
			
            if complete == false 				# reduce time step size concentration change criterion not satisfied
				dt *= knobs.dt_decrease
				assert(dt > knobs.dt_min)
			end
			
		end	
			
        # update values
        t += dt
		for (i, nd) in enumerate(node)
			nd.u += du[i]
		end
        
        # update time step
        dt *= knobs.dt_increase
        dt = min(dt, knobs.dt_max, t_end - t)	
	
	end
	
	# extract concentration from node type and return
	for (i, nd) in enumerate(node)
		u[i] = nd.u
	end
	return u[2:num_nodes-1]
	
end


# main script

function BoundedDiff(mode::Int64, t::Float64, num_pts::Int64)

	println("Reading parameter values.")
	solute, bar = GetParams()
	R = R_coeff(solute, bar)

	println("Delineating column.")
	dL = bar.L/num_pts
	x = zeros(Float64,num_pts)
	for i = 1:num_pts
		x[i] = dL * (i-0.5)
	end
	
	if mode==1
		println("Solving via analytical solution.")
		C = zeros(Float64,num_pts)
		D_eff = solute.D * (bar.phi	- bar.theta)
		for i = 1:num_pts
			C[i] = C_a(x[i], t, D_eff, R, bar)
		end		
	else
		println("Solving via numerical solution.")
		knobs = GetKnobs()										# read in numerical knobs
		C = C_n(x, t, solute.D, R, bar, knobs, num_pts) 		# solve by finite difference method
	end
	
	println("Writing output.")
	fname = "results.csv"
	csvfile = open(fname,"w")
	line_out = "x" * "," * "C"
	println(csvfile, line_out)	
	for i = 1:num_pts
		line_out = string(x[i]) * "," * string(C[i])	
		println(csvfile,line_out)
	end
	close(csvfile)		
	
	println("Done.")

end


### run script as BoundedDiff(mode, t, num_pts)

# mode==1 --> analytical solution; mode==2 --> numerical solution
# t = model time
# num_pts = number of sample locations (or nodes, if numerical) along column for output

BoundedDiff(2, 30.0, 100) 				# change these parameter values as warranted before running script