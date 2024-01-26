import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from matplotlib.ticker import FuncFormatter

### Functions regarding the setting of the parameters by inputs
def variable_builder():
    '''This function permits to set by inputs the variables necessary for the other functions, considering the two thermostats configuration.

    Parameters: ///

    Returns: the simulation interval of time, the bar lenght, the left thermostat temperature, 
    the right thermostat temperature, the temperature of the bar at the initial state.'''

    print('I need the following initial condition:')
    time_   = float(input("Simulation time (seconds) --> "))
    lenght_   = float(input("Lenght of the bar (metres): "))
    temperature_left_ = float(input("Temperature in the first extreme of the bar (Kelvin) --> "))
    temperature_right_ = float(input("Temperature in the second extreme of the bar (Kelvin) --> "))
    temperature_bar_  = float(input("Starting temperature of the other bar points (Kelvin) --> "))
    return time_,lenght_,temperature_left_,temperature_right_,temperature_bar_
def variable_builder_well():
    '''This function permits to set by inputs the variables necessary for the other functions, considering the three themostats configuration.

    Parameters: ///

    Returns: the simulation interval of time, the bar lenght, the left thermostat temperature, 
    the right thermostat temperature, the temperature of the bar at the initial state, the position and the temperature of the third thermostat.'''
    
    print('I need the following initial condition:')
    time_   = float(input("Simulation time (seconds) --> "))
    lenght_   = float(input("Lenght of the bar (metres) --> "))
    temperature_left_ = float(input("Temperature in the first extreme of the bar (Kelvin) --> "))
    temperature_right_ = float(input("Temperature in the second extreme of the bar (Kelvin) --> "))
    temperature_bar_  = float(input("Starting temperature of the other bar points (Kelvin) --> "))
    temperature_well_  = float(input("Temperature of the well (Kelvin) --> "))
    well_position_   =   float(input("Position of the well (metres) -->"))
    return time_,lenght_,temperature_left_,temperature_right_,temperature_bar_,temperature_well_,well_position_


### Functions to build the initial state of the bar
def bar_builder (temperature_left, temperature_right, temperature_bar, dim_X=100, dim_t=100):
    '''This function build the initial state of the simulation bar, assigning to the extremities fixed temperatures equal to the thermostat temperatures,
    and assigning to the rest of the space points the initial bar temperature.
    
    Parameters: left thermostat temperature(double),right thermostat temperature(double),bar initial temperature(double).
    
    Returns: the initial state of the bar temperature and the values of the two thermostats kept costant in time (array dim_X x dim_t).'''

    bar=np.zeros(shape=(dim_X, dim_t))
    for m in range(dim_X):
        bar[m][0]= temperature_bar
        bar[m][1]= temperature_bar
    for n in range(dim_t):
        bar[0][n]= temperature_right
        bar[-1][n]= temperature_left  
    return bar

def bar_builder_well (temperature_left, temperature_right, temperature_bar, well_position, temperature_well, lenght, dim_X=100, dim_t=100):
    '''This function build the initial state of the simulation bar, assigning to the extremities fixed temperatures equal to the thermostat temperatures,
    and assigning to the rest of the space points the initial bar temperature. It also set the position of the third thermostat in the bar and its fixed 
    temperature
    
    Parameters: left thermostat temperature(double),right thermostat temperature(double),bar initial temperature(double),the third thermostat position(double),
    the third thermostat temperature(double), the lenght of the bar(double)
    
    Returns: the initial state of the bar temperature and the values of the three thermostats kept costant in time (array dim_X x dim_t). '''

    bar=np.zeros(shape=(dim_X, dim_t))
    well_position=int(well_position/lenght*dim_X)
    for m in range(dim_X):
        bar[m][0]= temperature_bar
        bar[m][1]= temperature_bar
    for n in range(dim_t):
        bar[0][n]= temperature_right
        bar[-1][n]= temperature_left
        bar[well_position][n]= temperature_well
    return bar


### Functions applied to the initial bar state to obtain the evolution of the temperature bar profile in time, using a finite difference method
def DFF(lenght, time, bar, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the DuFortFrankel method, considering the two thermostat configuration.
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or of the bar_builder_well functions (array dim_X x dim_t).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    del_X          = lenght/dim_X
    del_t          = time/dim_t
    s = linear_diffusion*del_t/(pow(del_X,2))
    for n in range(dim_t-2):
        for m in range(dim_X-2):
            bar[m+1][n+2]= ((1-2*s)/(1+2*s)*bar[m+1][n]) + ((2*s)/(1+2*s)*(bar[m][n+1]+bar[m+2][n+1]))
    return bar

def DFF_well(lenght, time, bar, well_position, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the DuFortFrankel method, considering the three thermostat configuration.
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or bar_builder_well functions (array dim_X x dim_t)
    and the position of the third thermostat (double).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]    
    del_X          = lenght/dim_X 
    del_t          = time/dim_t
    s = linear_diffusion*del_t/(pow(del_X,2))
    well_position=int(well_position/lenght*dim_X)
    for n in range(dim_t-2):
        for m in range(dim_X-2):
            if m+1 != well_position:
                bar[m+1][n+2]= ((1-2*s)/(1+2*s)*bar[m+1][n]) + ((2*s)/(1+2*s)*(bar[m][n+1]+bar[m+2][n+1]))
    return bar

def C_N(lenght, time, bar, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the Crank-Nicolson method, considering the two thermostat configuration.
    
    Parameters: the bar lenght (double), the simulation interval of time(double), 
    the initial state of the bar equals to the return of the bar_builder or bar_builder_well functions (array dim_X x dim_t).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    del_X          = lenght/dim_X
    del_t          = time/dim_t
    matrix=np.zeros(shape=(dim_X,dim_X))
    s = linear_diffusion*del_t/(pow(del_X,2))
    m=0
    for l in range(dim_X):
        matrix[0][0]= 1
        matrix[-1][-1]=1
        if l != 0 and l != dim_X-1:
            matrix[l][m]= s+1
            matrix[l][m-1]= -s/2
            matrix[l][m+1]= -s/2
        m+=1
    solution= np.zeros(shape=(dim_X))
    for n in range(dim_t-1):
        for i in range(dim_X):
            solution[-1]=bar[-1][0]
            solution[0]=bar[0][0]
            if i != dim_X-1 and i != 0:
                solution[i]= (s/2) * (bar[i-1][n]+bar[i+1][n]) + (1-s)*bar[i][n]
        temperatura= np.linalg.solve(matrix, solution)
        for m in range(dim_X):
            if m < dim_X-1:
                bar[m+1][n+1]= temperatura[m+1]
    return bar

def C_N_well(lenght, time, bar, well_position, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the DuFortFrankel method, considering the three thermostat configuration. 
    It is really important taking in consideration the stability of this method. Considering the introduction of the third thermostat, 
    the C_N method decreases its stability and to mantain acceptable results it is necessary decrease the dimension of the temperature array 
    when the bar_builder_well function is called (usually set equal to dim_X_=100,dim_t_=100 as defoult).
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or bar_builder_well functions (array dim_X x dim_t)
    and the position of the third thermostat (double).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]   
    del_X          = lenght/dim_X
    del_t          = time/dim_t
    well_position=int(well_position/lenght*dim_X)
    matrix=np.zeros(shape=(dim_X,dim_X))
    m=0
    for l in range(dim_X):
        matrix[0][0]= 1
        matrix[well_position][well_position]=1
        matrix[-1][-1]=1
        if l != 0 and l != dim_X-1 and l != well_position:
            matrix[l][m]= pow(del_X,2)+linear_diffusion*del_t
            matrix[l][m-1]= -del_t*linear_diffusion/2
            matrix[l][m+1]= -del_t*linear_diffusion/2
        m+=1
    solution= np.zeros(shape=(dim_X))
    for n in range(dim_t-1):
        for i in range(dim_X):
            solution[0]=bar[0][0]
            solution[-1]=bar[-1][0]
            solution[well_position]=bar[well_position][0]
            if i != dim_X-1 and i != 0 and i != well_position:
                solution[i]=del_t*bar[i-1][n]+(pow(del_X,2)-2*del_t)*bar[i][n]+del_t*bar[i+1][n]
        temperatura= np.linalg.solve(matrix, solution)
        for m in range(dim_X):
            if m+1 < dim_X-1 and m+1 != well_position:
                bar[m+1][n+1]= temperatura[m+1]
    return bar

def R_K(lenght, time, bar, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the Runge-Kutta method, considering the two thermostat configuration.
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or of the bar_builder_well functions (array dim_X x dim_t).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    del_X          = lenght/dim_X 
    del_t          = time/dim_t
    s              = linear_diffusion*del_t/pow(del_X,2)
    A = s/2 
    B = 1-(3/2)*s
    C = 2*(s-1) 
    for n in range(dim_t-1):
        for m in range(dim_X):
            if m == 1:
                bar[m][n+1] = bar[m][n] + s*( bar[m-1][n]*B+bar[m][n]*C+bar[m+1][n]*B+bar[m+2][n]*A)
            if m > 1 and m < dim_X-2:
                bar[m][n+1] = bar[m][n] + s*(bar[m-2][n]*A+bar[m-1][n]*B+bar[m][n]*C+bar[m+1][n]*B+bar[m+2][n]*A)
            if m == dim_X-2:
                bar[m][n+1] = bar[m][n] + s*(bar[m-2][n]*A+bar[m-1][n]*B+bar[m][n]*C+bar[m+1][n]*B)
    return bar

def R_K_well(lenght, time, bar, well_position, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the DuFortFrankel method, considering the three thermostat configuration.
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or bar_builder_well functions (array dim_X x dim_t)
    and the position of the third thermostat (double).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    del_X          = lenght/dim_X 
    del_t          = time/dim_t
    s              = linear_diffusion*del_t/pow(del_X,2)
    well_position=int(well_position/lenght*dim_X)
    A = s/2 
    B = 1-(3/2)*s
    C = 2*(s-1) 
    for n in range(dim_t-1):
        for m in range(dim_X):
            if m != well_position:
                if m == 1:
                    bar[m][n+1] = bar[m][n] + s*( bar[m-1][n]*B+bar[m][n]*C+bar[m+1][n]*B+bar[m+2][n]*A)
                if m > 1 and m < dim_X-2:
                    bar[m][n+1] = bar[m][n] + s*(bar[m-2][n]*A+bar[m-1][n]*B+bar[m][n]*C+bar[m+1][n]*B+bar[m+2][n]*A)
                if m == dim_X-2:
                    bar[m][n+1] = bar[m][n] + s*(bar[m-2][n]*A+bar[m-1][n]*B+bar[m][n]*C+bar[m+1][n]*B)
    return bar

def FTCS(lenght, time, bar, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the Forward Time Centered Space method, considering the two thermostat configuration.
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or of the bar_builder_well functions (array dim_X x dim_t).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    del_X          = lenght/dim_X 
    del_t          = time/dim_t
    s              = linear_diffusion*del_t/pow(del_X,2)
    for n in range(dim_t-1):
        for m in range(dim_X-2):
            bar[m+1][n+1]= bar[m+1][n]+ s*(bar[m][n]-2*bar[m+1][n]+ bar[m+2][n])
    return bar

def FTCS_well(lenght, time, bar, well_position, linear_diffusion=0.00002):
    '''This function permits to resolve the problem by the DuFortFrankel method, considering the three thermostat configuration.
    
    Parameters: the bar lenght (double),the simulation interval of time (double), 
    the initial state of the bar equal to the return of the bar_builder or bar_builder_well functions (array dim_X x dim_t)
    and the position of the third thermostat (double).
    
    Returns: the evolution of the bar temperature profile during the simulation (array dim_X x dim_t).'''

    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    del_X          = lenght/dim_X
    del_t          = time/dim_t
    s              = linear_diffusion*del_t/pow(del_X,2)
    well_position=int(well_position/lenght*dim_X)
    for n in range(dim_t-1):
        for m in range(dim_X-2):
            if m+1 != well_position:
                bar[m+1][n+1]= bar[m+1][n]+ s*(bar[m][n]-2*bar[m+1][n]+ bar[m+2][n])
    return bar


### Functions to visualize the results of the precedent methods

def graph_time_t(lenght, time, istant, bar):
    '''This functions draws a plot of the bar temperature profile in a specific istant of the simulation.
    
    Parameters: lenght of the bar (double), time of the simulation (double), istant of the simulation to be represented (double), 
    the evolution of the bar temperature profile during the simulation equal to the return of the functions which apply a finite difference methods (array dim_X x dim_t)
    
    Returns:///'''

    istant=int((istant/time)*bar.shape[1])
    ax=plt.axes()
    X = np.linspace(0,lenght,bar.shape[0])
    bar_time_t=np.zeros(bar.shape[0])
    for i in range(bar.shape[0]):
        bar_time_t[i]=bar[i][istant]
    ax.plot(X,bar_time_t)
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel('Space (m)')
    ax.set_title('Temperature profile at time '+str(istant/bar.shape[1]*time)+' s')
    plt.show()

def plot_3D (lenght, time, bar):
    '''This functions draws a 3D plot of the bar temperature development along the bar during the entire simulation time.
    
    Parameters: lenght of the bar (double), time of the simulation(double), the evolution of the bar temperature profile during the simulation 
    equal to the return of the functions which apply a finite difference methods (array dim_X x dim_t).
    
    Returns:///'''

    plt.figure()
    dim_X=bar.shape[0]
    dim_t=bar.shape[1]
    ax = plt.axes(projection ="3d")
    gridx , gridy = np.meshgrid (range(dim_t),range(dim_X))

    ### Functions to format well the values on the axes of the plot
    def custom_scale_formatter_t(value, tick_number):
        '''The function is necessary to visualize the real values on the time axes of the graph, and not the indices of the array. It works just a rescaling.'''

        return value * time / dim_t
    def custom_scale_formatter_X(value, tick_number):
        '''The function is necessary to visualize the real values on the space axes of the graph, and not the indices of the array. It works just a rescaling.'''

        return value * lenght / dim_X
    
    ax.xaxis.set_major_formatter(FuncFormatter(custom_scale_formatter_t))
    ax.yaxis.set_major_formatter(FuncFormatter(custom_scale_formatter_X))
    ax.plot_wireframe (gridx, gridy, bar, cstride=2, rstride=2,linewidth=0.5,cmap='viridis')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Space (m)')
    ax.set_zlabel('Temperature (K)')
    ax.set_title('3D simulation reproduction')
    plt.show()

def plot_evolution(lenght, time, bar):
    '''By this function it is possible generating an animation that represents the temperature profile of the bar at each istant of the simulation.
    
    Parameters: lenght of the bar (double), time of the simulation (double), the evolution of the bar temperature profile during the simulation  
    equal to the return of the functions which apply a finite difference methods (array dim_X x dim_t)

    Returns: ///'''

    frame_selection = np.arange(0, bar.shape[1], 2)
    def funct_predict_t(istant):
        '''This function it is necessary to collect the different frames of the animation, at each istant. 
        A function like this it is expected by the ani.FuncAnimation.
        '''

        ax=plt.axes()
        X = np.linspace(0,lenght,bar.shape[0])
        bar_istant=np.zeros(bar.shape[0])
        for m in range(bar.shape[0]):
            bar_istant[m]=bar[m][istant]
        ax.plot(X,bar_istant)
        plt.text(lenght/3,np.amax(bar), 'Clock: '+ str(round(istant/100*time,2)) +' s', fontsize=12, color='black')
        ax.set_ylabel('Temperature (K)')
        ax.set_xlabel('Space (m)')
        ax.set_title('Animated Simulation')
        if istant==np.amax(frame_selection):
            plt.close()
        def on_close(event):
            '''
            This function permits to fig.canvas.mpl_connect to stop the programme running when the animation is interrupted.
            '''

            plt.close()
        fig.canvas.mpl_connect('close_event', on_close)
        return ax
    fig = plt.figure()
    anim = ani.FuncAnimation(fig,funct_predict_t,frames=frame_selection,interval=100,repeat=False)
    plt.show()


### Functions to visualize the comparison between the different finite difference methods

def methods_comparison(lenght, time, bar):
    ''' This function permits to compare the different methods results at the end of the simulation, by plotting the four temperature profiles in the same graph.
    It puts in evidence which approach is the most stable, for example. It takes in consideration the configuration with just 2 thermostats.
    
    Parameters: lenght of the bar (double), time of the simulation (double),the evolution of the bar temperature profile during the simulation 
    equal to the return of the functions which apply a finite difference methods (array dim_X x dim_t).

    Results: ///'''

    bar_DFF  = DFF(lenght,time,bar)
    bar_C_N  = C_N(lenght,time,bar)
    bar_R_K  = R_K(lenght,time,bar)
    bar_FTCS = FTCS(lenght,time,bar)
    ax=plt.axes()
    X = np.linspace(0,lenght,bar.shape[0])
    ax.plot(X,bar_DFF[:,-1],'co',label='DFF')
    ax.plot(X,bar_C_N[:,-1],'g^',label='C_N')
    ax.plot(X,bar_R_K[:,-1],'y--',label='R_K')
    ax.plot(X,bar_FTCS[:,-1],'b-',label='FTCS')
    plt.xlabel('Space (m)')
    plt.ylabel('Temperature (K)')
    plt.legend(loc='upper left')
    plt.show()

def methods_comparison_well(lenght, time, bar, well_position):
    ''' This function permits to compare the different methods results at the end of the simulation, by plotting the three temperature profiles in the same graph.
    It puts in evidence which approach is the most stable, for example. It takes in consideration the configuration with just 3 thermostats.
    
    Parameters: lenght of the bar (double), time of the simulation (double),the evolution of the bar temperature profile during the simulation 
    equal to the return of the functions which apply a finite difference methods (array dim_X x dim_t), and the position of the third thermostat (double).

    Results: ///'''

    bar_DFF  = DFF_well(lenght,time,bar,well_position)
    bar_R_K  = R_K_well(lenght,time_,bar,well_position)
    bar_FTCS = FTCS_well(lenght,time,bar,well_position)
    ax=plt.axes()
    X = np.linspace(0,lenght,bar.shape[0])
    ax.plot(X,bar_DFF[:,-1],'co',label='DFF')
    ax.plot(X,bar_R_K[:,-1],'y--',label='R_K')
    ax.plot(X,bar_FTCS[:,-1],'b-',label='FTCS')
    plt.xlabel('Space (m)')
    plt.ylabel('Temperature (K)')
    plt.legend(loc='upper left')
    plt.show()

