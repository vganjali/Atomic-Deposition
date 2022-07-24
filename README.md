## Basic Physics: L-J potential
The Lennard-Jones potential (L-J potential) is a mathematical model to approximate the interaction between a pair of neutral atoms or molecules.
$V_{LJ}=4\epsilon[(\frac{\sigma}{r})^{12}-(\frac{\sigma}{r})^{6}]=\epsilon[(\frac{r_{m}}{r})^{12}-2(\frac{r_{m}}{r})^{6}]$
where
- $\epsilon$: Depth of the potential well
- $\sigma$: Distance at wich $V_{LJ}=0$
- $r_m$: Distance at wich $V_{LJ}$ is minimum ($=-\epsilon$)
- ($r_m=2^{1/6}\sigma\approx1.122\sigma$)
- $r^{-12}$: The repulsive term, describes _Pauli repulsion_ at short ranges due to overlapping electrin orbitals
- $r^{-6}$ The attractive term, describes attraction at long ranges (for example, _Van de Waals force_)
- The same equation is used in _many-particle_ environment
## Algorithm Flowchart
![image](https://user-images.githubusercontent.com/3451891/180628036-f5ebe4eb-5005-4f4d-b9af-f44d51f2eb51.png)
## Key Considerations
### Bottom Layer
- A well-arranged bottom layer is assumed to be pre-formed when the simulation starts
- Exerts attractive/repulsive force on incoming particles according to L-J potential
- Provides a general model for the simulation
### Adaptive Time-step
- Particle slows down when L-J potential changes sharply
- Increases the stablity around the equilibrium point
- $r_{k+1}=r_{k}+\Delta t \nu_{k}+{\Delta t}^2\frac{F_{r_k}}{m}$
- $r_{k+1}=\nu_{k}+\frac{\Delta t}{2m}(F_{r_k}+F_{r_{k+1}})$
- $\Delta t_{k+1}={1}/{[\nu_k \times exp(r_c-r_n)]}$
- ![image](https://user-images.githubusercontent.com/3451891/180628689-55ee6a8d-0911-42f5-9350-632b8240069d.png)
- $r_c$ is critical distance to look for neighbors
- $r_n$ is the distance from the closest neighbor

### Timeout
- A counter is the triggered whenever the first change in the sign of the force $F=-\frac{\partial V_{LJ}}{\partial r}$ is detected
- After a reasonable number of time-steps, timeout happens and the particle reaches the equilibrium
- Simplifies complexity due to large force from the bottom surface
- ![image](https://user-images.githubusercontent.com/3451891/180628981-b93c798e-aefe-426c-98f0-afc8b6939bfc.png)

### Velocity Damping
- Velocity of the particle is dampled if there is a change in sign of velocity or force
- Velocity sign change indicates collision, and damping emulates inelastic collision
- Force sign change indicates crissing the equilibrium and damping velocity around it simplifies the motion
- ![image](https://user-images.githubusercontent.com/3451891/180628811-794cfc2f-8ed6-4389-ad1f-988defebfb52.png)

### Grids
- **1**/**0** indicate presence/absence of old atoms
- Submatrix dimension $(2 critical distance+1) \times (2 critical distance+1)$
- ![image](https://user-images.githubusercontent.com/3451891/180628870-8df63cb4-53af-4197-a47f-a6db53fb5fd6.png)

### Force Submatrix
- At the beginning, components of _F_ are calculated for each point within the submatrix
- A new particle uses a look-up table (LUT) of pre-calculated force values
- This increases the computation speed even when number of particles is large
- ![image](https://user-images.githubusercontent.com/3451891/180629030-543205f8-3060-4557-bb19-6ff8c130d890.png)


### Real Physical Parameters

- For _Cu_ atoms:
- ![image](https://user-images.githubusercontent.com/3451891/180628937-726dbce3-8b9c-4347-9367-a8a6a5204a98.png)


![image](https://user-images.githubusercontent.com/3451891/180629104-1f860d90-c5b9-4d6b-aa85-c049deb5f4c9.png)
![image](https://user-images.githubusercontent.com/3451891/180629108-0cfbdc5f-2bcf-40f1-9655-9e500cb1115c.png)
![image](https://user-images.githubusercontent.com/3451891/180629110-369719a0-ad57-47ab-aff1-663362ff229e.png)



