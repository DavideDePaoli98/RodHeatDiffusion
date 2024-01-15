# Temperature Profile Development Along an Iron Rod

<p align="center">
The diffusion equation, or heat equation, is a partial differential equation used to describe several phenomena: for example, in chemistry to model the diffusion of the species concentration, or in physics, to characterize the temperature development in a spatial domain, defined a specific interval of time and some initial conditions. Considering an iron tod in contact with a heat source, it is possible to represent the diffusion equation as follow: 
</p>

<p align="center">
<img src="images\q_equation_1_dimension.png" alt="Testo alternativo" width="220" height="80">
</p>

<p align="center">
To resolve the equation by a computational way, it is necessary discretizing the equation by the finite difference methods. Derivatives can be approximated with difference dependent expression, where the differences are between the values assumed by the function at two distinct but close domain points:

<p align="center">
<img src="images\derivate_approximation.png" alt="Testo alternativo" width="240" height="70">
</p>

<p align="center">
The programme proposed in this project is able to simulate the variation of temperature along an iron tod, in contact with two or three thermostats, during a given interval of time. The programme implements four finite difference methods, which can be applied on two different configuration. The first one considers two thermostats in contact with the iron tod extremities, meanwhile the second considers also another thermostat in contact with a third poit of the tod. Finite difference methods, used during the computing solution to the problem, are Crank-Nicolson, Runge-Kutta 2, and, Forward Time Centered Space and DuFort-Frankel.

## Finite differences methods
### Crank-Nicolson
<p align="center"> 
Crank-Nicolson is a second-order in time method and it is implicit. In particular, it is a combination of Forward Euler's and Back Euler's methods, and is represented by the scheme:

<p align="center">
<img src="images\C_N_eq.png" alt="Testo alternativo" width="550" height="70">
</p>
<p align="center">
where i indicates the space point and j the time istant. S can be found as:
<img src="images\S_eq.png" alt="Testo alternativo" width="120" height="30">

<p align="center">
In this case, it is necessary the resolution of a system of equations at the same time, following the matricial form:

<p align="center">
<img src="images\C_N_matrices.png" alt="Testo alternativo" width="800" height="150">
</p>

### Forward Time Centered Space
<p align="center">
Forward Time Centered Space is a method to solving differential equations of the first-order in time and explicit, which exploits the Forward Euler's method concerning the time domain and the central difference about the space. The scheme follows:
<p align="center">
<img src="images\F_T_C_S_eq.png" alt="Testo alternativo" width="320" height="60">
</p>

### DuFort-Frankel
<p align="center">
DuFort-Frankel method exploits the central difference in both space and time. It is a second-order method in time. His scheme, applied to the heat
equation, becomes:

<p align="center">
<img src="images\DF_F_eq.png" alt="Testo alternativo" width="360" height="110">
</p>

### Runge_kutta
<p align="center">
The Runge-Kutta method, chosen in this project, is also known as the Middle Point Method. It is a two-stage, explicit, second-order Runge-Kutta method, 
characterized by the following scheme:
</p>
<p align="center">
<img src="images\R_K_eq.png" alt="Testo alternativo" width="500" height="60">,<img src="ABC_R_K_eq.png" alt="Testo alternativo" width="250" height="50">.
</p>

### Stability problem:
<p align="center">
One of the most important limit of the finite difference methods is their stability. It is in general correleted to high values of **s**: if **s** reaches too high values, the results obtained by the programme are not realistic and they are characterized by large anomalous oscillationts.But each method has a specific range of stability and, in some case, it is enough large to permit a good stability with almost any initial real states. To give an idea, here it is reported a series of graphs which show the stability power of the Crank-Nicolson method respect the other:
</p>
<p align="center">
<img src="images\stability.png" alt="Testo alternativo" width="1000" height="400">
</p>
<p align="center">
A: s=0,32 , B: s=1,04 , C: s=2,48 , D: s=2,64 , E: s=8,00 , F: s=12,00
</p>

## The Programme 
<p align="center">
The programme permits to applied one of the four finite difference methods to an initial rod state, computing the evolution in time of the temperature rod profile. There steps necessary to obtain a simulation and its visualization are three:

* **1**: it is required chosing the configuration of the problem (2 or three thermostats) and fixing the initial state conditions. Particularly, the temperature of each thermostat, the lenght of the bar, the time interval of the simulation. It is possible changing the value of the linear diffusion coefficient D, considering a variation of the material, directly from the code (the iron value is set as default);

* **2**: then, one of the possible methods has to be applied to the bar to   obtain an array that represents the develop of the temperature in each position at each istant of time;

* **3**: the results can be represented by three graphic ways. 
    * At first, it is possible plotting the temperature profile at a specific istant of the simulation.
    * The second option is visualizing the simulation by an animation that shows the temperature profile changing in time, step by step, till the end of the time intervall given to the programme.
    * Third, it is also allowed to generate a plot 3D that permits to observe the shape of the entire process by a unique graph. 

Another possibility is to generate overlap istant plots generating from differt methods, usefull to compare their stability (which methods are more stable and which not).

\
The project is divided as follow:

[The DiffEqLibrary](miasottocartella/): this file collects all the function written to the programme. They are necessary: to build the initial configuration of the rod, to compute the temperature profile development starting from the initial condition and with the various methods, to build the graphic representatios of the results.\
[The Interface](interface.py): by running this code it is possible generating the simulations without working directly with the DiffEqLibrary functions. It is studied to guide people step by step, and it is not necessary programming to use it.\
[The Simulation example](simulation_example.py): here it is reported a simple example concerning the steps to use correctly the library proposed. To be sure that the library is installed correcly, it is suggested starting from these code to verify it works well.\
[The Test file](testing_file): in this file all the functions of the library are tested, working with different parameters ad control that each functions returns mantain the characteristics expected.

Under, some examples of 3D are reported, obtained by DuFort Frankel running the [The Simulation example](simulation_example.py):
</p>
<p align="center">
<img src="images\simulation_3D.png" alt="Testo alternativo" width="800" height="450">
</p>
<p align="center">
<img src="images\simulation_3D_well.png" alt="Testo alternativo" width="800" height="450">
</p>
