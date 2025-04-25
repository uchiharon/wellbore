#boundary Import necessary libraries
using Gridap                # Main Gridap package for finite element analysis
using Gridap.Geometry       # For mesh and geometry handling
using Gridap.FESpaces       # For finite element spaces
using Gridap.MultiField     # For coupled multi-physics problems
using Gridap.Io             # For input/output operations
using Gridap.Fields         # For field operations
using Gridap.TensorValues   # For tensor operations
using Gridap.ODEs           # For time-dependent problems
using Gridap.CellData       # For cell data operations and projection
using WriteVTK              # For VTK file output (visualization)
using GridapGmsh            # For Gmsh mesh integration

# ============================================================================
# PROBLEM DESCRIPTION
# ============================================================================
# This is a plane strain poroelasticity formulation
# In plane strain, we assume εzz = 0 (no strain in z-direction)
# but σzz ≠ 0 (stress in z-direction can exist)
# This is appropriate for modeling soil/rock layers where the z-dimension is 
# constrained but stresses can develop in that direction

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

# Material properties
E = 20.0e6         # Young's modulus (Pa)
nu = 0.2          # Poisson's ratio
B = 0.8           # Biot coefficient (coupling between fluid pressure and solid stress)
M = 1.0e9         # Biot modulus (Pa) - related to fluid and solid compressibility
k = 1.0e-3        # Permeability (m^2) - how easily fluid flows through the medium
mu = 1.0e-3       # Fluid viscosity (Pa·s)

# Loading conditions
Pb = 31.5e6
p0 = 20.0e6

# Time stepping parameters
T = 0.0005          # Final time (s)
num_steps = 3   # Number of time steps
dt = T / num_steps # Time step size (s)

# ============================================================================
# DERIVED MATERIAL PROPERTIES
# ============================================================================
# Calculate Lamé parameters for plane strain formulation
# For plane strain, we use the same Lamé parameters as in 3D
lambda = E * nu / ((1 + nu) * (1 - 2 * nu))  # First Lamé parameter (Pa)
mu = E / (2 * (1 + nu))                      # Second Lamé parameter (shear modulus) (Pa)
k_mu = k / mu                                # Hydraulic conductivity (permeability/viscosity)

dirichlet_tags = ["top_bottom", "wellbore"]
#####################################
######## ADD YOU CODE HERE ##########
#####################################
