# Pore-Network Modeling of Polymer Flows in Porous media

This project is focused on developing efficient and accurate computational methods for simulating the behavior of polymer fluids in complex porous structures. By incorporating advanced rheological models, we aim to capture the intricate flow dynamics and flow-particle interactions that are crucial in various geological and engineering applications.

Please contact Dr. Si Suo email: s.suo@imerial.ac.uk.

## Key Features

- **Pore Network Modeling**  
  Our code implements pore-scale simulations to study polymer flow through interconnected pores. The approach considers both the complex geometry of the porous medium and the nonlinear rheology of fluids.
  For network extraction, please refer to [https://github.com/ImperialCollegeLondon/porescale](https://github.com/ImperialCollegeLondon/porescale).

- **Non-Newtonian Behavior**  
  Unlike Newtonian fluids, which exhibit linear viscosity, polymer fluids generally exhibit the shear-thinning behavior. This repository accommodates a range of non-Newtonian constitutive models to better capture these complex fluid dynamics.
  
- **Flow-particle interaction force**  
  We also developed a bespoke module to calculate flow-particle interaction force based on the above PNM, which helps to understand how shear-thinning features impact the grain drags in granular materials.

- **Virtual porous media**  
  We use the discrete element method to generate sphere packings, surrogating a wide range of porous media. The provided package can generate a random sphere packing according to a given porosity, particle number, and particle size distribution (PSD).  

## Getting Started

### Prerequisites

To use the code in this repository, you will need:

- Matlab 2022+

## How to Cite
Si Suo, Sajjad Foroughi, Martin J. Blunt, Catherine O'Sullivan, (2025), Pore-Network Modeling of Polymer Flow in Porous Media. Computers and Geotechnics, 182, 107142.
DOI: [https://doi.org/10.1016/j.compgeo.2025.107142](https://doi.org/10.1016/j.compgeo.2025.107142)
