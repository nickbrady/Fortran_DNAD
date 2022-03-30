# Fortran_DNAD
Implement DNAD (Automatic Differentiation) into Newman's `BAND` subroutine

This code combines 3 items already developed in the literature:
1. John Newman's `BAND` subroutine (Appendix C, page 619 of Electrochemical Systems 3rd Ed. by John Newman and Karen E. Thomas-Alyea (2004))
2. `fillmat` found in Chapter 7 (page 129) of the thesis: [Experimental and numerical investigation of mass transfer in electrochemical systems by J. Deliang Yang (Columbia University, 1997)](https://clio.columbia.edu/catalog/1987854?counter=1)
3. And the dual number automatic differentiation Fortran module `dnadmod` (originally published as ["DNAD, a Simple Tool for Automatic
Differentiation of Fortran Codes Using Dual Numbers" by Wenbin Yu and Maxwell Blair](https://www.sciencedirect.com/science/article/pii/S0010465513000027?casa_token=MpXIh34txb0AAAAA:vf9mYSrbAU3VNKE9MYdLnQkd2OpTSa2AW0D5sN9FNbCI9fkhPZw-UXEcbR_4-CYoKAwEXgXmivA)); `dnadmod` was originally forked from [joddlehod](https://github.com/joddlehod/dnad)

The use of automatic differentiation significantly simplifies the process of linearizing systems of differential equations. In addition this is done without incurring a large performance loss and maintains numerical accuracy to within machine precision.

# Example Problems
Some example problems are included to help the user understand how models can be developed within this framework

General Equations <br>
Analytical Equation: <br>
* ∂cᵢ/∂t = -∇⋅𝐍ᵢ + Rᵢ    <br>

Finite Volume (Control Volume): <br>
* ΔV ∂cᵢ/∂t = (Aₓᵢ⋅𝐍ᵢ - Aₓₒ⋅𝐍ₒ) + ΔV ⋅ Rⱼ    <br>

***
## Simple Diffusion equation:
* ∂c/∂t = D∇²c
* 𝐍 = -D ∇c
* R = 0
* BC-WEST : cₒ = 1.0
* BC-EAST : cₒ = 0.0


`Flux_(1) = -diff * dc0dx`  <br>
`Rxn_(1) = 0.0`             <br>
`Accum_(1) = c0/delT`       <br>
`BC_WEST_(1) = c0 - 1.0`    <br>
`BC_EAST_(1) = c0 - 0.0`    <br>

***
## Diffusion-Reaction
* ∂c/∂t = D∇²c
* 𝐍 = -D ∇c
* R = -kᵣₓ ⋅ c
* BC-WEST : cₒ = cbulk
* BC-EAST : 𝐍ₒ = 0

`Flux_(1) = -diff * dc0dx`  <br>
`Rxn_(1) = k_Rxn * c0`      <br>
`Accum_(1) = c0/delT`       <br>
`BC_WEST_(1) = c0 - cbulk`  <br>
`BC_EAST_(1) = flux_temp(1) - 0.0`  <br>

***
## Battery Electrode Equations
* ϵ ∂cₒ/∂t = D∇²cₒ + a iᵣₓ / F
* (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F
* 0 = -∇⋅𝐢₁ - a iᵣₓ
* 0 = -∇⋅𝐢₂ + a iᵣₓ

<br>

* (1) ϵ ∂cₒ/∂t = D∇²cₒ + a iᵣₓ / F
* 𝐍ₒ = -ϵ * (D_Li * ∇cₒ + z_Li * u_Li * cₒ F ∇Φ₂)
* Rₒ =  a iᵣₓ / F
* BC-WEST : cₒ = cbulk
* BC-EAST : 𝐍ₒ = 0

`Flux_(1) = -porosity * (diff_Li * dc0dx + z_1 * u_1 * c0 * Fconst * dPhi_2dx)` <br>
`Rxn_(1) = +volumetric_surface_area * i_rxn / Fconst` <br>
`Accum_(1) = (porosity) * c0/delT` <br>
`BC_WEST_(1) = c0 - cbulk` <br>
`BC_EAST_(1) = flux_temp(1) - 0.0` <br>

* (2) (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F
* 𝐍ₓ = 0
* Rₓ = -a iᵣₓ / F
* BC-WEST : (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F
* BC-EAST : (1-ϵ) ∂cₓ/∂t = - a iᵣₓ / F

`Flux_(2) = 0.0` <br>
`Rxn_(2) = -volumetric_surface_area * i_rxn / Fconst` <br>
`Accum_(2) = (1.0 - porosity) * c_x/delT` <br>
`BC_WEST_(2) = accum_temp(2) - rxn_temp(2)` <br>
`BC_EAST_(2) = accum_temp(2) - rxn_temp(2)` <br>


* (3) 0 = -∇⋅𝐢₁ - a iᵣₓ
* 𝐢₁ = -(1-ϵ) σ ∇Φ₁
* Rₓ =  a iᵣₓ
* BC-WEST : 𝐢₁ = 0                (Solid-State Current = 0)
* BC-EAST : Φ₂ = 0               (arbitrary ref Voltage)

`Flux_(3) = -(1 - porosity) * sigma * dPhi_1dx` <br>
`Rxn_(3) = -volumetric_surface_area * i_rxn` <br>
`Accum_(3) = 0.0` <br>
`BC_WEST_(3) = flux_temp(3) - 0.0` <br>
`BC_EAST_(3) = flux_temp(3) - applied_current_A` <br>


* (4) 0 = -∇⋅𝐢₂ + a iᵣₓ
* 𝐢₂/F = ∑ᵢ (zᵢ 𝐍ᵢ)
* Rₓ = -a iᵣₓ
* BC-WEST : 𝐢₁ = i_applied     (All current carried in solid state)
* BC-EAST : 𝐢₂ = 0             (Solution Current = 0)

`Flux_(4) = -porosity * ( Fconst * (z_1 * diff_Li + z_2 * diff_PF6) * dc0dx + (z_1**2 * u_1 + z_2**2 * u_2) * Fconst**2 * c0 * dPhi_2dx )` <br>
`Rxn_(4) = +volumetric_surface_area * i_rxn` <br>
`Accum_(4) = 0.0` <br>
`BC_WEST_(4) = Phi_2 - 0.0` <br>
`BC_EAST_(4) = flux_temp(4) - 0.0` <br>

------------------------------------------------------------------------------------------------------------------------
## Cooling Fluid in Pipe
* ∂(ρc T)/∂t = -∇⋅(ρcT⋅𝐯) + h (T - Tₐ)
* ∂ρ/∂t = -∇⋅(ρ𝐯)
* BC-WEST :
  * T = Tᵢ
  * 𝐯 = vᵢ
* BC-EAST :
  * ∇T = 0      (Temperature does not change after exiting)
  * ∇⋅(ρ𝐯) = 0  (ρ∇𝐯 + 𝐯⋅∇ρ = 0 ; ∇ρ = ∂ρ/∂T ∇T)

(1) Energy balance <br>
`Flux_(1) = c_heat_cap * density_ * Temp * vel` <br>
`Rxn_(1) = -h_heat_transfer * (Temp - T_ambient)` <br>
`Accum_(1) = c_heat_cap * density_ * Temp/delT` <br>
`BC_WEST_(1) = Temp - T_in` <br>
`BC_EAST_(1) = dTdx - 0.0` <br>

(2) Continuity equation <br>
`Flux_(2) = density_ * vel` <br>
`Rxn_(2) = 0.0` <br>
`Accum_(2) = density_/delT` <br>
`BC_WEST_(2) = vel - vel_in` <br>
`BC_EAST_(2) = density_%dx(1)*dTdx * vel + density_*dveldx - 0.0` <br>


# Code Structure
Because Fortran necessitates that modules appear before subsequent modules that depend on them, `include` statements have been utilized to maintain this rigid ordering, while keeping the core subroutines and modules hidden from the general user.
(Previously, the python script `SortFortranModules.py` was written to re-order the modules and then run the programs.) <br>
If one is not interested in using a python script to run the Fortran programs, one can elect to run the `.f95` files from the terminal by first compiling the program: `gfortran -fdefault-real-8 -O3 _file_name_` and then executing the program: `./a.out`

### Core_Subs_Mods
This folder contains all the subroutines and modules that make up the numerical engine to solve coupled non-linear partial-differential-equations. These subroutines and modules are written generally (so they can be used to solve a wide variety of problems) and therefore rarely need to be modified. Segregating these subroutines and modules from the other code enables us to more easily maintain "master" versions of the essential code and eliminates clutter for the general user, who can now more easily focus on the code specific to their particular system.
