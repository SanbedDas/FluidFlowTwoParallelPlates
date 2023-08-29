# Fluid Flow Between Two Parallel Plates
**Fluid Flow Simulation using the Semi-Implicit Method for Pressure Linked Equations (SIMPLE)**

The **Semi-Implicit Method for Pressure Linked Equations (SIMPLE)** is a widely utilized numerical technique for solving the Navier-Stokes equations that govern incompressible fluid flow. This method holds immense significance in various engineering applications due to its ability to accurately simulate fluid flow patterns, heat distribution, and pressure within such systems.

**Key Features:**
- The algorithm employs a staggered grid system wherein pressure and velocity are located at distinct grid points, ensuring stability and accuracy in simulations.
- The procedure consists of initialization, momentum prediction, pressure correction, and velocity correction steps, each contributing to an accurate representation of fluid behavior.
- Stability is achieved by implicitly treating momentum equations while explicitly handling pressure equations, allowing for more reliable simulations and enabling selection of appropriate time steps.

**Implementation Details:**
- The code focuses on solving the Navier-Stokes equations for an incompressible fluid within a 2D rectangular domain.
- The domain is discretized using a uniform grid with nx and ny grid points in the x and y directions.
- Velocity components (u and v) in the x and y directions, pressure (P), fluid density (ρ), and dynamic viscosity (µ) constitute the governing parameters.

**Abstract:**
The **SIMPLE algorithm** plays a pivotal role in simulating incompressible fluid flow, heat distribution, and pressure variations. Its staggered grid framework, involving pressure and velocity at distinct locations, enhances simulation stability. The initialization, momentum prediction, pressure correction, and velocity correction steps contribute to accurate fluid behavior representation. Noteworthy is the approach's implicit treatment of momentum equations and explicit handling of pressure equations, striking a balance between stability and accuracy. The method's application is demonstrated through solving the Navier-Stokes equations within a 2D domain using a uniform grid discretization. This README introduces the SIMPLE algorithm's significance, implementation details, and its capacity to effectively simulate incompressible fluid flow phenomena.
