
# Integrators for Newtonian motion

module Integrator

using MyVec

export PlutoSim, Verlet, RungeKutta, SIA4


#const G = 6.67398e-11 # G in basic SI units
const G = 1.48811382e-34 # G in units of AU and day

function potential{T<:Real}(mass::Vector{T}, r::Vector{IVector{T}})
    N = size(r,1)
    out = zeros(IVector{T}, N)
    for i = 1:N
        for j = i+1:N
            d = r[j]-r[i]
            F = G / norm(d)^3
            a = d * F
            out[i] += mass[j] * a
            out[j] -= mass[i] * a
        end
    end
    out
end

# Magic constants for the SIA4 algorithm
const CONST_B1 =  0.5153528374311229364
const CONST_B2 = -0.085782019412973646
const CONST_B3 =  0.4415830236164665242
const CONST_B4 =  0.1288461583653841854
const CONST_C1 =  0.1344961992774310892
const CONST_C2 = -0.2248198030794208058
const CONST_C3 =  0.7563200005156682911
const CONST_C4 =  0.3340036032863214255

# 4th order symplectic
function SIA4{T<:Real}(mass::Vector{T}, r0::Vector{IVector{T}}, v0::Vector{IVector{T}}, h::Real)    
    v1 = v0 + CONST_B1 * h * potential(mass, r0)
    r1 = r0 + CONST_C1 * h * v1
    
    v2 = v1 + CONST_B2 * h * potential(mass, r1)
    r2 = r1 + CONST_C2 * h * v2
    
    v3 = v2 + CONST_B3 * h * potential(mass, r2)
    r3 = r2 + CONST_C3 * h * v3
    
    v4 = v3 + CONST_B4 * h * potential(mass, r3)
    r4 = r3 + CONST_C4 * h * v4
    return (r4, v4)
end



PlutoSim(N::Integer, h::Real) = PlutoSim(N,h,100,Verlet)
PlutoSim(N::Integer, h::Real, M::Integer) = PlutoSim(N,h,M,Verlet)


function PlutoSim(N::Integer, h::Real, M::Integer, f::Function)
    r = Array(IVector{Float64}, 6)
    v = Array(IVector{Float64}, 6)
#    mass = zeros(6)
    
    # Set masses
    const mass = [1.988435e30, 1.8988e27, 5.685e26, 8.6625e25, 1.0278e26, 1.314e22]
    
    ## Set starting points (AU) and velocity vectors (AU/day)
    # Sun
    r[1] = [-9.001119860001391E-04, -2.518569331216286E-03, -5.035993149103694E-05]
    v[1] = [ 6.249336869691155E-06, -7.202874592454979E-07, -1.368159714329672E-07]
    
    # Jupiter
    r[2] = [ 9.637969147947649E-01,  4.988694611226984E+00, -4.236695676623525E-02]
    v[2] = [-7.500214649417347E-03,  1.791097950812253E-03,  1.603945618163088E-04]
    
    # Saturn
    r[3] = [-7.901119375041307E+00, -5.800240922607197E+00,  4.152972681975652E-01]
    v[3] = [ 2.998916003562968E-03, -4.511136746229775E-03, -4.061079792754435E-05]

    # Uranus
    r[4] = [ 1.985459827926987E+01,  2.801201477916861E+00, -2.468196041542967E-01]
    v[4] = [-5.782111510118402E-04,  3.711193104939768E-03,  2.130487704096845E-05]

    # Neptune
    r[5] = [ 2.664809431198752E+01, -1.375183029235338E+01, -3.309310983483438E-01]
    v[5] = [ 1.419153776960717E-03,  2.808496228357725E-03, -9.056396274034148E-05]
    
    # Pluto 
    r[6] = [5.303572815981638E+00, -3.190544944470667E+01, 1.879962140526394E+00]
    v[6] = [3.156551097775996E-03, -1.132909607403660E-04, -9.009397396467777E-04]
    
    Nsaves = fld(N,M)
    result = zeros(Nsaves+1, 36)
    println(size(result))
    println("Begin.")
    t0 = time()
    t_start = t0
    for i = 1:Nsaves
        # Print time statistics
        t1 = time()
        tot = t1 - t_start
        delta = t1-t0
        est = (tot / (i-1)) * (Nsaves - i)
        @printf("%5d / %5d: last = %7.3f, total = %7.3f, remain = %7.3f\n", i, Nsaves, delta, tot, est)
        t0 = t1
        # Save results
        for p = 1:6            
            result[i, 6*(p-1)+1] = r[p,1]
            result[i, 6*(p-1)+2] = r[p,2]
            result[i, 6*(p-1)+3] = r[p,3]
            result[i, 6*(p-1)+4] = v[p,1]
            result[i, 6*(p-1)+5] = v[p,2]
            result[i, 6*(p-1)+6] = v[p,3] 
        end
        # Actual computation
        for j = 1:M
            r,v = SIA4(mass, r, v, h)
        end
    end
    # Save last results
    for p = 1:6            
        result[Nsaves+1, 6*(p-1)+1] = r[p,1]
        result[Nsaves+1, 6*(p-1)+2] = r[p,2]
        result[Nsaves+1, 6*(p-1)+3] = r[p,3]
        result[Nsaves+1, 6*(p-1)+4] = v[p,1]
        result[Nsaves+1, 6*(p-1)+5] = v[p,2]
        result[Nsaves+1, 6*(p-1)+6] = v[p,3] 
    end
    
    println("Writing output...")
    writecsv("output.txt", result)
#    println("Running plot script...")
#    run(`python plotorbits.py`)
end


end # module