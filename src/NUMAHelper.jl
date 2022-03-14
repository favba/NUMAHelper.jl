module NUMAHelper

export read_cons_file, read_last_error, read_dx_resolution, read_error, read_dt

"""
Reads a .cons file returning a named tuple of vectors (itime,time,mass,energy)
"""
function read_cons_file(filename::AbstractString)

  open(filename) do f
    N = parse(Int,readline(f)) + 1

    itime = Vector{Int}(undef,N)
    time = Vector{Float64}(undef,N)
    mass = Vector{Float64}(undef,N)
    energy = Vector{Float64}(undef,N)

    for i=1:N
      s = split(readline(f))
      itime[i] = parse(Int,s[1])
      setindex!.((time,mass,energy),parse.(Float64,(s[2:4]...,)),i)
    end
    
    return (itime=itime, time=time, mass=mass, energy=energy)
  end

end

"""
Reads the dx resolution used in the simulation
"""
function read_dx_resolution(infile)
	s = readlines(infile)
	i = findfirst(x->(occursin("Grid Resolution: dx",x)),s)
	parse(Float64,split(s[i])[end])
end

"""
Reads the time step used in the simulation
"""
function read_dt(infile)
	s = readlines(infile)
	i = findfirst(x->(occursin("dt time_initial time_final time_restart time_scale =",x)),s)
	parse(Float64,split(s[i])[7])
end

"""
Reads the error of the solution at time "time"
"""
function read_error(infile::AbstractString,time::Real)
  s = reverse(readlines(infile))
  error_names = (:l1,:l2,:l8,:rms)
  i = findfirst(s) do line
	    if occursin("itime time dt =",line)
		    return parse(Float64,split(line)[6]) ≈ time ? true : false
	    else
		    return false
	    end
    end
  i += 8
  j = i
  err = (;zip(error_names,split(s[j])[7:10] .|> x->parse(Float64,x))...)

  names_func = function (l)
	  j = i-l+1
  	d = split(s[j])
	  Symbol(d[1][1:end-1])
  end

  variable_names = ntuple(names_func,5)

  errors_func = function (l)
	  j = i-l+1
    err = (;zip(error_names,split(s[j])[7:10] .|> x->parse(Float64,x))...)
  end

  errors = ntuple(errors_func,5)
  r = (;zip(variable_names,errors)...)
end

"""
Reads the last error outputed in the simulation
"""
function read_last_error(instring::AbstractString)
  s = reverse(readlines(instring))
  error_names = (:l1,:l2,:l8,:rms)
  i = findfirst(x->occursin("Rho:    L1 L2 L8 RMS ",x),s)
  j = i
  err = (;zip(error_names,split(s[j])[7:10] .|> x->parse(Float64,x))...)

  names_func = function (l)
	  j = i-l+1
  	d = split(s[j])
	  return Symbol(d[1][1:end-1])
  end

  variable_names = ntuple(names_func,5)

  errors_func = function (l)
	  j = i-l+1
    return (;zip(error_names,split(s[j])[7:10] .|> x->parse(Float64,x))...)
  end

  errors = ntuple(errors_func,5)
  r = (;zip(variable_names,errors)...)
end

end # module
