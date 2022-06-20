function cart_to_geo(x,y,z)
    deg = 180/π
    longitude = atan(y,x)*deg # Between -180 and 180
    H2 = x^2 + y^2
    latitude = atan(z/(sqrt(H2)))*deg # Between -90 and 90
    height = sqrt(H2 + z^2)
    return longitude, latitude, height
end

function cart_to_geo!(lon,lat,h,x,y,z)
    @inbounds for i in eachindex(x)
        lon[i],lat[i],h[i] = cart_to_geo(x[i],y[i],z[i])
    end
end

function vec_cart_to_geo(u,v,w,x,y,z)
    phi = atan(y,x)
    theta = atan(sqrt(x^2 + y^2)/z)
  
    vz = v*cos(phi) - u*sin(phi)
    vm = u*cos(phi)*cos(theta) + v*sin(phi)*cos(theta) - w*sin(theta)
    vr = u*cos(phi)*sin(theta) + v*sin(phi)*sin(theta) + w*cos(theta) 
  
    return vz, vm, vr
end

function vec_cart_to_geo(u,v,w,longitute,latitude)
    rad = π/180
    phi = longitute*rad
    theta = (90-latitude)*rad
  
    vz = v*cos(phi) - u*sin(phi)
    vm = u*cos(phi)*cos(theta) + v*sin(phi)*cos(theta) - w*sin(theta)
    vr = u*cos(phi)*sin(theta) + v*sin(phi)*sin(theta) + w*cos(theta) 
  
    return vz, vm, vr
end

function vec_cart_to_geo!(vz::AbstractArray,vm::AbstractArray,vr::AbstractArray,u::AbstractArray,v::AbstractArray,w::AbstractArray,point...)
    @inbounds for i in eachindex(vz)
        vz[i],vm[i],vr[i] = vec_cart_to_geo(u[i],v[i],w[i],getindex.(point,(i,))...)
    end
end

function read_var!(output::AbstractArray,ncvar,I::AbstractArray)
    @inbounds for (ii,i) in enumerate(I)
        output[ii] = ncvar[i]
    end
    return output
end

function read_var(ncvar,I::AbstractArray)
    output = zeros(eltype(ncvar),length(I))
    return read_var!(output,ncvar,I)
end

function read_var(filename::AbstractString,var::AbstractString,I::AbstractArray)
    NetCDF.open(filename) do f
        read_var(f[var],I)
    end
end

function read_var!(output::AbstractArray,filename::AbstractString,var::AbstractString,I::AbstractArray)
    NetCDF.open(filename) do f
        read_var!(output,f[var],I)
    end
end

read_var(filename::AbstractString,var::AbstractString) = ncread(filename,var)
read_var!(output::AbstractArray,filename::AbstractString,var::AbstractString) = ncread!(filename,var,output)

function unique_z(z)
    set = SortedSet{Float32}()
    @inbounds for z_i in z
        Float32(z_i) ∉ set && push!(set,Float32(z_i))
    end
    return set
end

function unique_xy(x,y)
    set = SortedSet{NTuple{2,Float32}}()
    @inbounds for i in eachindex(x)
        xy = Float32.((x[i],y[i]))
        xy ∉ set && push!(set,xy)
    end
    return set
end

function return_index(set::SortedSet,v)
    i = 1
    for val in set
        val == v && break
        i+=1
    end
    return i
end

function Base.getindex(set::SortedSet,i::Integer)
    j = 0
    r = zero(eltype(set))
    @inbounds for val in set
        j+=1
        if j == i
            r = val
            break
        end
    end
    return r
end

function read_var_by_z_levels(filename::AbstractString,var::AbstractString,z)

    z_unique = unique_z(z)
    ny = length(z_unique)
    nx = count(x->(Float32(x) == z_unique[1]),z)
    result = zeros(nx,ny)
    
    @inbounds for (j,val) in enumerate(z_unique)
        ind = findall(x->(Float32(x) == val),z)
        read_var!(@view(result[:,j]),filename,var,ind)
    end

    return result
end

function compute_quantity_at_pressure_level(quantity::AbstractArray, plevel::Number, pressure::AbstractArray)
    result = zeros(size(pressure,1))
    @inbounds for i in axes(pressure,1)
        ind = -1
        for j in axes(pressure,2)
            if pressure[i,j] <= plevel
                ind = j
                break
            end
        end
        dp = pressure[i,ind] - pressure[i,ind-1]
        result[i] = ((pressure[i,ind]-plevel)/dp)*quantity[i,ind-1] + ((plevel-pressure[i,ind-1])/dp)*quantity[i,ind]
    end
    return result
end