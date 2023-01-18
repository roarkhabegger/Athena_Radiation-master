<comment>
problem   = Simulation with Turbulence, CII+Ly Cooling function
reference = 
configure =  python configure.py --prob=chem_H2_MHD --chemistry=H2_MHD -fft -mpi -hdf5 -b --cvode_path=/global/homes/k/kyuen2/cvode/cvode/ --lib_path=/global/homes/k/kyuen2/cvode/cvode/lib64/ --nscalars=2 --ccmd CC

<job>
problem_id = MP  # problem ID: basename of output filenames

<output1>
file_type  = hst         # History data dump
dt         = 0.001       # time increment between outputs

<output2>
file_type  = hdf5        # Binary data dump
variable   = prim        # variables to be output
dt         = 1.0         # time increment between outputs

<output3>
file_type  = rst         # History data dump
dt         = 5.0        # time increment between outputs

<time>
cfl_number = 0.3        # The Courant, Friedrichs, & Lewy (CFL) Number
nlim       = 1 # 10000000   # cycle limit
tlim       = 5.0 # 100.0      # time limit
integrator  = vl2       # time integration algorithm
xorder      = 2         # order of spatial reconstruction
ncycle_out  = 1         # interval for stdout summary info

<mesh>
nx1        = 32         # Number of zones in X1-direction
x1min      = -100.0     # minimum value of X1
x1max      = 100.0      # maximum value of X1 (unit: pc)
ix1_bc     = periodic   # inner-X1 boundary flag
ox1_bc     = periodic   # outer-X1 boundary flag

nx2        = 32         # Number of zones in X2-direction
x2min      = -100.0     # minimum value of X2
x2max      = 100.0      # maximum value of X2 (unit: pc)
ix2_bc     = periodic   # inner-X2 boundary flag
ox2_bc     = periodic   # outer-X2 boundary flag

nx3        = 32         # Number of zones in X3-direction
x3min      = -100.0     # minimum value of X3
x3max      = 100.0      # maximum value of X3 (unit: pc)
ix3_bc     = periodic   # inner-X3 boundary flag
ox3_bc     = periodic   # outer-X3 boundary flag

<meshblock>
nx1        = 16
nx2        = 16
nx3        = 16

<hydro>
gamma           = 1.666666666667 # gamma = C_p/C_v
iso_sound_speed = 1.57363272           # equavalent to sqrt(gamma*p/d) for p=0.1, d=1 (From examples) (unit: km s^-1)

<cr>
vmax             =2.99792458e5
sigmaPerp = 2.99792458e1
sigmaParl = 2.99792458e1
src_flag = 1
vs_flag = 1

<problem>
turb_flag  = 3        # 2  )1 for decaying, 2 (impulsive) or 3 (continuous) for driven turbulence
nH         = 3.0      # The density of each pixel (Unit: mH cm^-3)
invbeta    = 0.2      # inverse of Plasma beta (B^2/(2P))
angle      = 0.0      # The B-field angle to y-direction
invbetaCR  = 0.2      # inverse of CR beta (P_C/P)
vx         = 0.0      # The rms of velocity

cooling    = 1
s_init_H   = 1       # H+ aboundance per H
s_init_H2  = 0       # H+ aboundance per H

nH_c       = 3.0
b0_c       = 1.72382827 
pres_c     = 7.42895978
crpres_c   = 1.48579196
rad_c      = 20.0
wid_c      = 0.1

# Note in CR Athena implementation turbulence parameters fall in problem block
dedt       = 1e7      # 200 )Energy injection rate (for driven) or Total energy (for decaying)
nlow       = 0        # cut-off wavenumber at low-k
nhigh      = 12       # cut-off wavenumber at high-k
expo       = 2.0      # power-law exponent
tcorr      = 0.5      # correlation time for OU process (both impulsive and continuous)
dtdrive    = 0.2      # time interval between perturbation (impulsive)
rseed      = -1       # if non-negative, seed will be set by hand (slow PS generation)

<chemistry>
#chemistry solver parameters, please ignore it right now
reltol     = 1.0e-2    #relative tolerance, default 1.0e-2
abstol     = 1.0e-12   #absolute tolerance, default 1.0e-12
#abstol_H2 = 1.0e-10   #set absolute tolerance to a certain species
user_jac   = 0         #flag for whether use user provided Jacobian. default false/0
maxsteps   = 200000    #maximum number of steps in one integration. default 10000
h_init     = 1e-4      #first step of first zone. Default 0/CVODE algorithm.
output_zone_sec = 0     #output diagnostic
#code units
unit_density_in_nH = 1.
unit_length_in_cm = 3.086e+18
unit_vel_in_cms = 1.0e5

