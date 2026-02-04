# 1D Heat Diffusion Equations for Urban Surfaces

This document describes the mathematical formulation for 1D heat diffusion in urban surfaces, including pervious and impervious roads, roofs, and walls.

## Pervious Road Heat Diffusion

### Governing Equation

The 1D heat diffusion equation for pervious road soil layers is:

$$\frac{\partial}{\partial t}\left[C_v(z) T(z,t)\right] = 
  \frac{\partial}{\partial z}\left[k(z) \frac{\partial T}{\partial z}\right] + S(z,t)$$

where:
- $T(z,t)$ is temperature at depth $z$ and time $t$ [K]
- $C_v(z) = c_{solid}(1-\theta_{sat})\Delta z + (h_{ice} c_{pice} + h_{liq} c_{pliq})$ 
  is volumetric heat capacity [J/K]
- $k(z)$ is thermal conductivity [W/m/K]
- $S(z,t)$ is the source term (only at surface) [W/m²]

### Spatial Discretization

Using finite volumes with layer centers at $z_{c,j}$ and interfaces at $z_{i,j+1/2}$:

$$(C_v)_j \Delta z_j \frac{dT_j}{dt} = 
  k_{j+1/2}\frac{T_{j+1} - T_j}{\Delta z_{j+1/2}} - 
  k_{j-1/2}\frac{T_j - T_{j-1}}{\Delta z_{j-1/2}} + S_j$$

where $\Delta z_{j+1/2} = z_{c,j+1} - z_{c,j}$ is the distance between layer centers.

### Surface Boundary Condition (Top Layer, $j=0$)

The surface energy flux enters through the source term at the new time level:

$$S_0^{n+1} = h^{n+1} - \frac{\partial h}{\partial T}^{n}\left( T_0^{n+1} - T_0^n\right)$$

where the net surface energy flux is:

$$h = \overrightarrow{SW}_{net} - \overrightarrow{LW}_{net} - H_g - \lambda_v(E_{soil} + E_{tran})$$

with:
- $\overrightarrow{SW}_{net}$ = net shortwave radiation absorbed [W/m²]
- $\overrightarrow{LW}_{net}$ = net longwave radiation emitted [W/m²]  
- $H_g$ = sensible heat flux to atmosphere [W/m²]
- $E_{soil}$ = evaporation from soil [kg/m²/s]
- $E_{tran}$ = transpiration from vegetation [kg/m²/s]
- $\lambda_v$ = latent heat of vaporization [J/kg]

The derivative term for linearization:

$$\frac{\partial h}{\partial T} = \frac{\partial \overrightarrow{SW}_{net}}{\partial T} - \frac{\partial \overrightarrow{LW}_{net}}{\partial T} - \frac{\partial H_g}{\partial T} - \lambda_v\left(\frac{\partial E_{soil}}{\partial T} + \frac{\partial E_{tran}}{\partial T}\right)$$

where:
- $\frac{\partial \overrightarrow{SW}_{net}}{\partial T} = 0$ (shortwave radiation is independent of surface temperature)
- $\frac{\partial \overrightarrow{LW}_{net}}{\partial T} = 4\varepsilon\sigma T^3$ (longwave radiation derivative from Stefan-Boltzmann law)
- $\frac{\partial H_g}{\partial T}$, $\frac{\partial E_{soil}}{\partial T}$, and $\frac{\partial E_{tran}}{\partial T}$ are described in the surface flux documentation

The effective surface depth uses a tuning factor $\alpha = 0.34$:

$$\Delta z_{eff,0} = \frac{1}{2}(\Delta z_{1,0} + \alpha \Delta z_{2,0})$$

### Bottom Boundary Condition (Layer $j = N-1$)

Zero heat flux at the bottom:

$$k_{N-1/2}\frac{\partial T}{\partial z}\bigg|_{z=z_{bottom}} = 0$$

### Thermal Properties (Pervious Road)

#### Thermal Conductivity (Johansen 1975 Model)

Layer thermal conductivity:

$$k_j = K_e \cdot k_{sat} + (1 - K_e) \cdot k_{dry}$$

where the Kersten number $K_e$ depends on soil state:

$$K_e = \begin{cases}
  \max(0, \log_{10}(S_w) + 1) & T \geq T_{frz} \text{ (unfrozen)} \\
  S_w & T < T_{frz} \text{ (frozen)}
\end{cases}$$

with degree of saturation:

$$S_w = \frac{\theta_{liq} + \theta_{ice}}{\theta_{sat}}$$

The saturated thermal conductivity (geometric mean):

$$k_{sat} = k_{minerals} \cdot k_{water}^{f_l \theta_{sat}} \cdot k_{ice}^{(1-f_l)\theta_{sat}}$$

where $f_l = \theta_{liq}/(\theta_{liq} + \theta_{ice})$ is the liquid fraction.

#### Interface Thermal Conductivity (Harmonic Mean)

$$k_{j+1/2} = \frac{k_j \cdot k_{j+1} \cdot (z_{c,j+1} - z_{c,j})}
  {k_j(z_{c,j+1} - z_{i,j+1/2}) + k_{j+1}(z_{i,j+1/2} - z_{c,j})}$$

### Temporal Discretization (Crank-Nicolson)

Using Crank-Nicolson with weight $\theta_{CN} = 0.5$ (trapezoidal rule):

$$\frac{T_j^{n+1} - T_j^n}{\Delta t} = 
  \theta_{CN}\left[\frac{F_{j+1/2}^{n+1} - F_{j-1/2}^{n+1}}{(C_v)_j \Delta z_j}\right] + 
  (1-\theta_{CN})\left[\frac{F_{j+1/2}^n - F_{j-1/2}^n}{(C_v)_j \Delta z_j}\right]$$

where $F_{j+1/2} = k_{j+1/2}(T_{j+1} - T_j)/\Delta z_{j+1/2}$ is the heat flux.

### Tridiagonal System

This yields the linear system: $a_j T_{j-1}^{n+1} + b_j T_j^{n+1} + c_j T_{j+1}^{n+1} = r_j$

#### Top Layer ($j=0$):

$$r_0 = T_0^n + \beta_0\left(h - \frac{\partial h}{\partial T}T_0^n 
  + \theta_{CN}F_{1/2}^n\right)$$

$$b_0 = 1 + (1-\theta_{CN})\beta_0\frac{k_{1/2}}{\Delta z_{1/2}} 
  - \beta_0\frac{\partial h}{\partial T}$$

$$c_0 = -(1-\theta_{CN})\beta_0\frac{k_{1/2}}{\Delta z_{1/2}}$$

where $\beta_j = \Delta t/((C_v)_j \Delta z_j)$ for top layer with $\Delta z_{eff,0}$.

#### Interior Layers ($0 < j < N-1$):

$$r_j = T_j^n + \beta_j \theta_{CN}(F_{j+1/2}^n - F_{j-1/2}^n)$$

$$a_j = -(1-\theta_{CN})\beta_j\frac{k_{j-1/2}}{\Delta z_{j-1/2}}$$

$$b_j = 1 + (1-\theta_{CN})\beta_j\left(\frac{k_{j+1/2}}{\Delta z_{j+1/2}} 
  + \frac{k_{j-1/2}}{\Delta z_{j-1/2}}\right)$$

$$c_j = -(1-\theta_{CN})\beta_j\frac{k_{j+1/2}}{\Delta z_{j+1/2}}$$

#### Bottom Layer ($j=N-1$):

$$r_{N-1} = T_{N-1}^n - \theta_{CN}\beta_{N-1}F_{N-1/2}^n + \beta_{N-1}F_{N-1/2}$$

$$a_{N-1} = -(1-\theta_{CN})\beta_{N-1}\frac{k_{N-3/2}}{\Delta z_{N-3/2}}$$

$$b_{N-1} = 1 + (1-\theta_{CN})\beta_{N-1}\frac{k_{N-3/2}}{\Delta z_{N-3/2}}$$

The system is solved using the Thomas algorithm (tridiagonal matrix algorithm).

## Impervious Road Heat Diffusion

*To be documented*

## Roof Heat Diffusion

*To be documented*

## Wall Heat Diffusion

*To be documented*

## References

- Johansen, O. (1975). Thermal conductivity of soils. PhD thesis, University of Trondheim, Norway.
- ELM Documentation: SoilTemperatureMod.F90
