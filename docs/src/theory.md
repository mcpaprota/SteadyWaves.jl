# Theoretical background

Here, a brief descrition of the plane wave problem and its solution based on a Fourier Approximation Method is presented.
We consider steady, periodic waves propagatin in water of constant depth. The fluid is inviscid and incompressible
and its motion is irrotational. Thus, The velocity vector is expressed by the scalar velocity potential function as $\mathbf{v}=\nabla\varPhi(x,z,t)$. The only restoring force is gravity, since surface tension is neglected.
A schematic view of the problem is presented in the figure below.



Under these assumptions, the fluid mass is conserved according to the Laplace equation

$$\frac{\partial^2\varPhi}{\partial x^2}+\frac{\partial^2\varPhi}{\partial z^2} = 0, \quad 0 \leq z \leq \eta,$$

while the Euler equations of motion are presented as the Cauchy-Lagrange intergral

$$\frac{\partial\varPhi}{\partial t} + \frac{1}{2}\left(\left(\frac{\partial\varPhi}{\partial x}\right)^2 + \left(\frac{\partial\varPhi}{\partial z}\right)^2\right) + gz + \frac{p}{\rho} = 0, \quad 0 \leq z \leq \eta,$$

known as the unsteady Bernoulli equation. In the equations $g$ is the acceleration due to gravity, while $\rho$ and $p$ are the fluid density and pressure, respectively.