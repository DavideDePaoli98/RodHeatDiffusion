# Temperature Profile Development Along an Iron Rod

The programme proposed in this project is able to simulate the variation of temperature along an iron tod during a given interval of time. The tod can be considered in contact with two or three thermostats (on the extremities and on another tod point). The programme implements four finite difference methods to compute the development of the temperature tod profile:  Crank-Nicolson, Runge-Kutta 2, Forward Time Centered Space and DuFort-Frankel. 

### How install the programme

To install the programme clone the repository [DiffEqLibrary](https://github.com/DavideDePaoli98/RodHeatDiffusion) and use pip:
```
git clone https://github.com/DavideDePaoli98/RodHeatDiffusion
cd RodHeatDiffusion
pip install -r requirements.txt
```
To have an idea about how using the library concerning the finite difference methods, an example is given in [the simulation file](simulation_example.py). 

### How using the programme 

Going into detail, the steps necessary to obtain a simulation and its visualization are three:
* **1**: first, it is necessary setting in the [setting parameters file](setting_parameters.txt) the values of the different initial state conditions. Particularly, the temperature of each thermostat, the lenght of the bar, the time interval of the simulation. It is possible changing the value of the linear diffusion coefficient D, considering a variation of the material, directly from the code (the iron value is set as default parameter). If you want to consider the case with just two thermostats, you can ignore the values well_position and temperature_well.

* **2**: secondly, using the functions of the [DiffEqLibrary](DiffEqLibrary/), it is possible building the initial tod state and applying one of the finite difference method on the tod, to obtain the resulting temperature profile at the end of the simulation. To set the parameters it is required using the function *parameters_setting*, to obtain the initial bar state the function *bar_builder* or *bar_builder_well*, to develop the temperature profile during the simulation the functions relative to the finite difference methods (for example *FTCS* and *FTCS_well* applied on the initial bar state).

* **3**: finally, the result can be represented by three graphic ways: 
    * a 2D plot (space and temperature) of the temperature profile at a specific istant of the simulation, using the function *graph_time_t*;
    * an 2D animation (space and temperature) showing the temperature profile changing in time till the end of the simulation, using *plot_evolution*;
    * a 3D plot (space, time, temperature) that permits to observe the shape of the entire process by a unique graph, using *plot_3D*. 

Another possibility is to generate an istant plots overlap comparing the differt methods, by the function *methods_comparison* or *methods_comparison_well*. 

</p>

## Finite differences methods applyied on the diffusion equation

The diffusion equation, or heat equation, is a partial differential equation used to describe several phenomena: for example, in chemistry to model the diffusion of the species concentration, or in physics, to characterize the temperature development in a spatial domain, defined a specific interval of time and some initial conditions. Considering an iron tod in contact with a heat source, it is possible to represent the diffusion equation as follow: 


<p align="center">
<img src="images\q_equation_1_dimension.png" alt="Testo alternativo" width="220" height="80">
</p>


where D is equal to the termic conduttivity over the specific volume of the material considered and its density. In the programme, D corresponds to the material iron set as default, but it is possible changing it directly from the code.
To resolve the equation by a computational way, it is necessary discretizing the equation by the finite difference methods. Derivatives can be approximated with difference dependent expression, where the differences are between the values assumed by the function at two distinct but close domain points:

<p align="center">
<img src="images\derivate_approximation.png" alt="Testo alternativo" width="240" height="70">
</p>


### Crank-Nicolson

Crank-Nicolson is a second-order in time method and it is implicit. In particular, it is a combination of Forward Euler's and Back Euler's methods, and is represented by the scheme:

<p align="center">
<img src="images\C_N_eq.png" alt="Testo alternativo" width="550" height="70">
</p>

where i indicates the space point and j the time istant. S can be found as:
<p align="center">
<img src="images\S_eq.png" alt="Testo alternativo" width="120" height="30">
</p>

In this case, it is necessary the resolution of a system of equations at the same time, following the matricial form:

<p align="center">
<img src="images\C_N_matrices.png" alt="Testo alternativo" width="800" height="150">
</p>

### Forward Time Centered Space

Forward Time Centered Space is a method to solving differential equations of the first-order in time and explicit, which exploits the Forward Euler's method concerning the time domain and the central difference about the space. The scheme follows:
<p align="center">
<img src="images\F_T_C_S_eq.png" alt="Testo alternativo" width="320" height="60">
</p>

### DuFort-Frankel

DuFort-Frankel method exploits the central difference in both space and time. It is a second-order method in time. His scheme, applied to the heat
equation, becomes:

<p align="center">
<img src="images\DF_F_eq.png" alt="Testo alternativo" width="360" height="110">
</p>

### Runge_kutta

The Runge-Kutta method, chosen in this project, is also known as the Middle Point Method. It is a two-stage, explicit, second-order Runge-Kutta method, 
characterized by the following scheme:

<p align="center">
<img src="images\R_K_eq.png" alt="Testo alternativo" width="500" height="60"> 
</p>
<p align="center">
<img src="images\ABC_R_K_eq.png" alt="Testo alternativo" width="250" height="50">
</p>

### Stability problem:

One of the most important limit of the finite difference methods is their stability. It is in general correleted to high values of s: if s reaches too high values, the results obtained by the programme are not realistic and they are characterized by large anomalous oscillationts.But each method has a specific range of stability and, in some case, it is enough large to permit a good stability with almost any initial real states. To give an idea, here it is reported a series of graphs which show the stability power of the Crank-Nicolson method respect the other:

<p align="center">
<img src="images\stability.png" alt="Testo alternativo" width="1000" height="400">
</p>
<p align="center">
A: s=0,32 , B: s=1,04 , C: s=2,48 , D: s=2,64 , E: s=8,00 , F: s=12,00
</p>

## The Directory organization and some example

The project is divided as follow:

[The DiffEqLibrary](DiffEqLibrary/): this file collects all the function written to the programme. They are necessary: to build the initial configuration of the rod, to compute the temperature profile development starting from the initial condition and with the various methods, to build the graphic representatios of the results.\
[The Simulation example](simulation_example.py): here it is reported a simple example concerning the steps to use correctly the library proposed. To be sure that the library is installed correcly, it is suggested starting from these code to verify it works well.\
[The Test file](testing.py): in this file all the functions of the library are tested, verifying that each functions return mantains the expected characteristics and controlling some specific cases.\
[The Requirements file](requirements.txt): the file collected the libraries necessary to the programme.\
[The Setting Parameters File](setting_parameters.txt): this file permits to set the initial condition of the simulation without operating directly on the code.

Under, some examples of 3D are reported, obtained by DuFort Frankel method:
</p>
<p align="center">
<img src="images\simulation_3D.png" alt="Testo alternativo" width="600" height="450">
</p>
<p align="center">
<img src="images\simulation_3D_well.png" alt="Testo alternativo" width="600" height="450">
</p>
