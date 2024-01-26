from DiffEqLibrary import DiffEqLibrary_ as fp

''' First, it is necessary giving to the programme,by the setting_parameters.txt file, the variables to build the initial state of the bar: 
the time of the simulation, the lenght of the bar, the three temperatures of the two thermostats and of the bar at the beginning of the simulation.'''
parameters=fp.parameters_setting()


''' Secondly, the funcion bar-build generates the initial configuration of the array representing the bar.''' 
initial_bar = fp.bar_builder(temperature_left,temperature_right,temperature_bar)

''' Third, the functions correlated to the finite difference methods are applied on the initial bar to obtain the evolution of the temperature along the bar in time.
In this case, the method DuFortFrankel is selected to compute the develop of the simulation.'''
final_bar   = fp.DFF(lenght,time,initial_bar)

# Then, it is possible to work on the representation of the simulation. 

''' The function graph_time_t permits to visualize the temperature profile in one istant of the simulation;
in this case, the istant chosen correspond to the entire time of the simulation over 2. '''
fp.graph_time_t(lenght,time,time/2,final_bar)

'''The function plot_3d represents the entire evolution of the bar temperature profile during the simulation time,
by a reproduction of a 3D plot (Time vs Space vs Temperature).'''
fp.plot_3D (lenght,time,final_bar)

''' The function plot_evolution represents the entire evolution of the bar temperature profile during the simulation time by an animation;
the clock permits to visualize the istant of each profile. '''
fp.plot_evolution(lenght,time,final_bar)
