# Pore-Network Modeling of Polymer Flows in Porous media

This project is focused on developing efficient and accurate computational methods for simulating the behavior of polymer fluids in complex porous structures. By incorporating advanced rheological models, we aim to capture the intricate flow dynamics that are crucial in various geological and engineering applications.

Please contact Dr. Si Suo email: s.suo@imerial.ac.uk.

## Key Features

- **Pore Network Modeling**  
  Our code implements pore-scale simulations to study polymer flow through interconnected pores. The approach considers both the complex geometry of the porous medium and the nonlinear rheology of fluids.
  For network extraction, please refer to https://github.com/ForoughiSajjad/pnextract.

- **Non-Newtonian Behavior**  
  Unlike Newtonian fluids, which exhibit linear viscosity, polymer fluids generally exhibit the shear-thinning behavior. This repository accommodates a range of non-Newtonian constitutive models to better capture these complex fluid dynamics.

- **Virtual porous media**  
  We use the discrete element method to generate sphere packings, surrogating a wide range of porous media. The provided package can generate a random sphere packing according to a given porosity, particle number, and particle size distribution (PSD).  

## Getting Started

### Prerequisites

To use the code in this repository, you will need:

- Matlab 2022+
