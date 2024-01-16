from DiffEqLibrary import DiffEqLibrary_ as fp
import numpy as np

def main():
    print(
     '\nHello, I am a simulation software called Joy. \nI am able to simulate the temperature profile of a iron bar in contact with thermostatos.\n'
    +'Particularly, I can resolve two different configurations:\n   --> two thermostas in contact with the extremities of the bar;\n '
    +'  --> three thermostats, two in contact with the extremities of the bar and one in another poit, chosen by you.\n' 
    +'\nEach situation can be resolved by four different approaches, or better four different finite-difference methods: \n'
    +'  a) Crank-Nicolson         b) Forward Time Centered Space\n'+'  c) DuFortFrankel          d) Runge-Kutta\n'
    +'\nThe results I can produce are divided in three categories:\n'
    +'  1) The temperature along the bar in a specific istant of time after the beginning of the simulation.\n'
    +'  2) The evolution of the temperature along the bar during an interval of time, starting at the beginning of the simulation.\n'
    +'  3) The 3D plot that represents the evolution of the temperature bar in time (Time,Space,Temperature).\n'
    +'\nTo process these results, I need some initial condition:\n'
    +'  The bar lenght, the simulation interval of time, the extremity thermostat temperatures,\n'
    +'  the initial bar temperature, the position and the temperature of the possible third thermostat.\n'
    +'\nThe default linear diffusion coefficient D is set considering an iron bar,\n'
    +'but it is possible changing it directly operating on the code.\n\n\n ')

    answer = int(input('First of all, do you want to study a two thermostats configuration(1) or a three thermostats configuration(2)?\n'))
    variables=[]
    if answer == 1:
        variables=fp.variable_builder()
        bar=fp.bar_builder(temperature_left_=variables[2],temperature_right_=variables[3],temperature_bar_=variables[4])
        bar=without_well(initial_bar_=bar,variables_=variables)
    elif answer == 2:
        variables=fp.variable_builder_well()
        bar=fp.bar_builder_well(temperature_left_=variables[2],temperature_right_=variables[3],temperature_bar_=variables[4],well_position_=variables[6],temperature_well_=variables[5],lenght_=variables[1])
        bar=with_well(initial_bar_=bar,variables_=variables)
    else:   
        while answer != 1 and answer !=2:
            answer = int(input('\nI do not understand your answer, you have to select the first option, digiting 1 and push Enter\n'
            +'or the second one,digiting 2 and push Enter.\n'+'Which option do you chose?\n'))
        if answer == 1:
            variables=fp.variable_builder()
            bar=fp.bar_builder(temperature_left_=variables[2],temperature_right_=variables[3],temperature_bar_=variables[4])
            bar=without_well(initial_bar_=bar,variables_=variables)
        elif answer == 2:
            variables=fp.variable_builder_well() 
            bar=fp.bar_builder_well(temperature_left_=variables[2],temperature_right_=variables[3],temperature_bar_=variables[4],well_position_=variables[6],temperature_well_=variables[5],lenght_=variables[1])
            bar=with_well(initial_bar_=bar,variables_=variables)
    if bar == 'quit':
        answer = 'quit'
    while answer != 'quit':
        print('How do you want to visualize the results?\n')
        answer=str(input('\nIf you want visualize the temperature profile of the bar in a specific istant, digit istant.\n'
        +'If you want to visualize the evolution of the temperature profile in time during the interval indicated by you, digit evolution.\n'
        +'If you want to visualize the 3D plot that represents the evolution of the temperature bar in time (Time,Space,Temperature), digit 3D.\n'
        +'If you want to quit, digit quit.'
        +'What is your choice?\n'))
        if answer == 'istant':
            istant = float(input('Which is the istant you want to represents?'))
            fp.graph_time_t(lenght_=variables[1],time_=variables[0],istant_=istant,bar_=bar)
        elif answer == 'evolution':
            fp.plot_evolution(lenght_=variables[1],time_=variables[0],bar_=bar)
        elif answer == '3D':
            fp.plot_3D(lenght_=variables[1],time_=variables[0],bar_=bar)
        else:   
            while answer != 'istant' and answer != 'evolution' and answer != '3D' and answer != 'quit':
                answer = str(input('\nI do not understand your answer, you have to select the first option, digiting istant and pushing Enter\n'
                +'or the second one,digiting evolution and pushing Enter, or the third, digiting 3D and pushing Enter\n'
                +'If you want to quit, digit quit.\n'
                +'Which option do you chose?\n'))
                if answer == 'istant':
                    istant = float(input('Which is the istant you want to represents?'))
                    fp.graph_time_t(lenght_=variables[1],time_=variables[0],istant_=istant,bar_=bar)
                elif answer == 'evolution':
                    fp.plot_evolution(lenght_=variables[1],time_=variables[0],bar_=bar)
                elif answer == '3D':
                    fp.plot_3D(lenght_=variables[1],time_=variables[0],bar_=bar)
def without_well(initial_bar_,variables_):
    answer = str(input( '\nWhich finite-differences method I have to use?\n a --> Crank-Nicolson\n'
        +'  b --> Forward Time Centered Space\n  c --> DuFortFrankel\n  d --> Runge-Kutta\n'
        +' Instead, if you prefer to observe the comparison between the different possible methods, digit comparison\n'))
    if answer == 'a':
        final_bar = fp.C_N(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
    elif answer == 'b':
        final_bar = fp.FTCS(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
    elif answer == 'c':
        final_bar = fp.DFF(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
    elif answer == 'd':
        final_bar = fp.R_K(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
    elif answer == 'comparison':
        fp.methods_comparison(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
        final_bar = 'quit' 
    else:   
        while answer != 'a' and answer != 'b' and answer != 'c' and answer != 'd' and answer != 'comparison':
            answer = str(input('\nI do not understand your answer, you have to select one of the four option,\n'
            +'digiting a and push Enter for the first, digiting b and push Enter for the second,\n'
            +'digiting c and push Enter for the third, or digiting d and push Enter for the fourth.\n'+'Which option do you chose?\n'))
        if answer == 'a':
            final_bar = fp.C_N(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
        elif answer == 'b':
            final_bar = fp.FTCS(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
        elif answer == 'c':
            final_bar = fp.DFF(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
        elif answer == 'd':
            final_bar = fp.R_K(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
        elif answer == 'comparison':
            fp.methods_comparison(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_)
            final_bar = 'quit' 
    return final_bar
def with_well(initial_bar_,variables_):
    answer = str(input( '\nWhich finite-differences method I have to use?\n'
        +'  b --> Forward Time Centered Space\n  c --> DuFortFrankel\n  d --> Runge-Kutta\n'
        +'\nInstead, if you prefer to observe the comparison between the different possible methods, digit comparison\n'
        +'(I suggest not using the Cranck-Nicolson method with this configuration, because there is a problem given by the default code setting.\n'))
    if answer == 'b':
        final_bar = fp.FTCS_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
    elif answer == 'c':
        final_bar = fp.DFF_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
    elif answer == 'd':
        final_bar = fp.R_K_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
    elif answer == 'comparison':
        fp.methods_comparison_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
        final_bar = 'quit'    
    else:   
        while answer != 'a' and answer != 'b' and answer != 'c' and answer != 'd' and answer != 'comparison':
            answer = str(input('\nI do not understand your answer, you have to select one of the four option,\n'
            +'digiting a and push Enter for the first, digiting b and push Enter for the second,\n'
            +'digiting c and push Enter for the third, or digiting d and push Enter for the fourth.\n'+'Which option do you chose?\n'))
        if answer == 'a':
            final_bar = fp.C_N_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
        elif answer == 'b':
            final_bar = fp.FTCS_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
        elif answer == 'c':
            final_bar = fp.DFF_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
        elif answer == 'd':
            final_bar = fp.R_K_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
        elif answer == 'comparison':
            fp.methods_comparison_well(lenght_=variables_[1],time_=variables_[0],bar_=initial_bar_,well_position_=variables_[6])
            final_bar = 'quit' 
    return final_bar

main()