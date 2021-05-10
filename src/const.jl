const pi = 3.14159265358979323846264338327   #Pi
const grav = 9.8                               #Gravitational acceleration (m / s^2)
const cp = 1004.0                             #Specific heat of dry air at constant pressure
const cv = 717.0                              #Specific heat of dry air at constant volume
const rd = 287.0                              #Dry air constant for equation of state (P=rho*rd*T)
const p0 = 1.e5                              #Standard pressure at the surface in Pascals
const C0 = 27.5629410929725921310572974482   #Constant to translate potential temperature into pressure (P=C0*(rho*theta)**gamma)
const gamma = 1.40027894002789400278940027894   #gamma=cp/Rd 
#Define domain and stability-related constants
const xlen = 2.e4    #Length of the domain in the x-direction (meters)
const zlen = 1.e4    #Length of the domain in the z-direction (meters)
const hv_beta = 0.25     #How strong to diffuse the solution: hv_beta \in [0:1]
const cfl = 1.50    #"Courant, Friedrichs, Lewy" number (for numerical stability)
const max_speed = 450        #Assumed maximum wave speed during the simulation (speed of sound + speed of wind) (meter / sec)
const hs = 2          #"Halo" size: number of cells beyond the MPI tasks's domain needed for a full "stencil" of information for reconstruction
const sten_size = 4          #Size of the stencil used for interpolation

#Parameters for indexing and flags
const NUM_VARS = 4           #Number of fluid state variables
const ID_DENS = 1           #index for density ("rho")
const ID_UMOM = 2           #index for momentum in the x-direction ("rho * u")
const ID_WMOM = 3           #index for momentum in the z-direction ("rho * w")
const ID_RHOT = 4           #index for density * potential temperature ("rho * theta")
const DIR_X = 1              #Integer constant to express that this operation is in the x-direction
const DIR_Z = 2              #Integer constant to express that this operation is in the z-direction

