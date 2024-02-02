from DiffEqLibrary import DiffEqLibrary_ as fp
import numpy as np
from hypothesis import given
import hypothesis.strategies as st


#Test to control if the different functions work well with general initial conditions

@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000))
def test_bar_builder(temperature_left,temperature_right,temperature_bar):
    
    '''Test to control if the return of the bar builder have fixed temperatures in the extremetes and
    along the initial bar (considering all the simulation time).'''
    
    bar=fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    assert all(i == temperature_right for i in bar[0,:])
    assert all(i == temperature_left for i in bar[-1,:])
    assert all(i == temperature_bar for i in bar[1:-1,0]) 
@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),temperature_well=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5))
def test_bar_builder_well(temperature_left,temperature_right,temperature_bar,temperature_well,lenght):
    
    '''Test to control if the return of the bar builder have fixed temperatures in the extremetes,
    in the well and along the initial bar (considering all the simulation time).'''
    
    well_position=lenght/2
    bar=fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    well_position=int(well_position/lenght*bar.shape[0])
    assert all(i == temperature_right for i in bar[0,:])
    assert all(i == temperature_left for i in bar[-1,:])
    assert all(i == temperature_well for i in bar[well_position,:])
    assert all(i == temperature_bar for i in bar[1:well_position,0])
    assert all(i == temperature_bar for i in bar[well_position+1:-1,0])


@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300))
def test_DFF(temperature_left,temperature_right,temperature_bar,lenght,time):

    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes 
    during the entire simulation time.'''
    
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.DFF(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_bar for i in final_bar[1:-1,0])
@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300),temperature_well=st.floats(min_value=0.01,max_value=100000))
def test_DFF_well(temperature_left,temperature_right,temperature_bar,lenght,time,temperature_well):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes
    and in the well during the entire simulation time. '''
    
    well_position=lenght/2
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.DFF_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])
    assert all(i == temperature_bar for i in final_bar[1:well_position,0])
    assert all(i == temperature_bar for i in final_bar[well_position+1:-1,0])


@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300))
def test_C_N(temperature_left,temperature_right,temperature_bar,lenght,time):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes 
    during the entire simulation time. '''
    
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.C_N(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300),temperature_well=st.floats(min_value=0.01,max_value=100000))
def test_C_N_well(temperature_left,temperature_right,temperature_bar,lenght,time,temperature_well):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes
    and in the well during the entire simulation time.'''
    
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.C_N_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])


@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300))
def test_R_K(temperature_left,temperature_right,temperature_bar,lenght,time):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes 
    during the entire simulation time.'''
    
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.R_K(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300),temperature_well=st.floats(min_value=0.01,max_value=100000))
def test_R_K_well(temperature_left,temperature_right,temperature_bar,lenght,time,temperature_well):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes
    and in the well during the entire simulation time. '''
    
    well_position=int(well_position/lenght*initial_bar.shape[0])
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.R_K_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])


@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300))
def test_FTCS(temperature_left,temperature_right,temperature_bar,lenght,time):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes 
    during the entire simulation time.'''
    
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.FTCS(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
@given(temperature_left=st.floats(min_value=0.01,max_value=100000),temperature_right=st.floats(min_value=0.01,max_value=100000),temperature_bar=st.floats(min_value=0.01,max_value=100000),lenght=st.floats(min_value=0.1,max_value=5),time=st.floats(min_value=1,max_value=300),temperature_well=st.floats(min_value=0.01,max_value=100000))
def test_FTCS_well(temperature_left,temperature_right,temperature_bar,lenght,time,temperature_well):
    
    '''Test to control if the return of the finite difference method functions have fixed temperatures in the extremetes
    and in the well during the entire simulation time. '''
    
    well_position=int(well_position/lenght*initial_bar.shape[0])
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.FTCS_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])


 
#Test to verify the programme results in specific cases

@given(temperature_left=st.just(300.0),temperature_right=st.just(300.0),temperature_bar=st.just(300.0),lenght=st.just(1.0),time=st.just(10.0))
def test_case1(temperature_left,temperature_right,temperature_bar,lenght,time):
    
    '''If the left thermostat temperature , the right thermostat temperature and the initial bar temperature are equals, 
    it is expected that the final temperature profile is equal to the initial one and that the temperature along the bar 
    is always equal to the temperature_left'''
    
    bar = fp.bar_builder(temperature_left,temperature_right,temperature_bar)
    final_bar = fp.DFF(lenght,time,np.copy(bar))
    assert np.all(final_bar == temperature_left)
    assert np.array_equal(bar,final_bar)
@given(temperature_left=st.just(300.0),temperature_right=st.just(300.0),temperature_bar=st.just(300.0),lenght=st.just(1.0),time=st.just(10.0),well_position=st.just(0.5),temperature_well=st.just(300.0))
def test_case1_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    
    '''If the the left thermostat temperature , the right thermostat temperature, the initial bar temperature 
    and the well temperature are equals, it is expected that the final temperature profile is equal to the initial 
    one and that the temperature along the bar is always equal to the temperature_left'''
    
    bar_well = fp.bar_builder_well(temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar_well = fp.DFF_well(lenght,time,np.copy(bar_well),well_position)
    assert np.all(final_bar_well == temperature_left)
    assert np.array_equal(bar_well,final_bar_well)

@given(temperature_left=st.just(300.0),temperature_right=st.just(300.0),temperature_bar=st.just(200.0),lenght=st.just(1.0),time=st.just(10.0))
def test_case2(temperature_left,temperature_right,temperature_bar,lenght,time):
    
    '''If the left  and the right thermostat temperatures are higher than the bar temperature, 
    it is expected that the final temperature profile values are always lower or equal than the thermostat 
    temperatures.'''
    
    bar = fp.bar_builder(temperature_left,temperature_right,temperature_bar)
    final_bar = fp.DFF(lenght,time,np.copy(bar))
    assert np.all(final_bar <= temperature_left)
    assert np.all(final_bar <= temperature_right)
@given(temperature_left=st.just(300.0),temperature_right=st.just(300.0),temperature_bar=st.just(200.0),lenght=st.just(1.0),time=st.just(10.0),well_position=st.just(0.5),temperature_well=st.just(170.0))
def test_case2_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    
    '''If the left  and the right thermostat temperatures are higher than the bar temperature and the well temperature, 
    it is expected that the final temperature profile values are always lower or equal than the thermostat 
    temperatures and higher than the well temperature.'''
    
    bar_well = fp.bar_builder_well(temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar_well = fp.DFF_well(lenght,time,np.copy(bar_well),well_position)
    assert np.all(final_bar_well <= temperature_left)
    assert np.all(final_bar_well <= temperature_right)
    assert np.all(final_bar_well >= temperature_well)

@given(temperature_left=st.just(300.0),temperature_right=st.just(200.0),temperature_bar=st.just(250.0),lenght=st.just(1.0),time=st.just(10.0))
def test_case3(temperature_left,temperature_right,temperature_bar,lenght,time):
    
    '''If the left thermostat temperature is higher than the bar temperature, and if the right thermostat 
    temperature is lower than the bar temperature, the temperature profile will decreases going from the left 
    to the right.'''
    
    bar = fp.bar_builder(temperature_left,temperature_right,temperature_bar)
    final_bar = fp.DFF(lenght,time,np.copy(bar))
    assert np.all(final_bar <= temperature_left)
    assert np.all(final_bar >= temperature_right)
    assert all(i <= i+1 for i in final_bar[:,-1])
@given(temperature_left=st.just(300.0),temperature_right=st.just(200.0),temperature_bar=st.just(250.0),lenght=st.just(1.0),time=st.just(10.0),well_position=st.just(0.5),temperature_well=st.just(230.0))
def test_case3_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    
    '''If the left thermostat temperature is higher than the bar temperature and the well temperature, 
    and if the right thermostat temperature is lower than the bar temperature and the well temperature, 
    the temperature profile will decreases going from the left to the right.'''
    
    bar_well = fp.bar_builder_well(temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar_well = fp.DFF_well(lenght,time,np.copy(bar_well),well_position)
    assert np.all(final_bar_well <= temperature_left)
    assert np.all(final_bar_well >= temperature_right)
    assert all(i <= i+1 for i in final_bar_well[:,-1])

@given(temperature_left=st.just(200.0),temperature_right=st.just(300.0),temperature_bar=st.just(300.0),lenght=st.just(1.0),time=st.just(10.0),well_position=st.just(0.5),temperature_well=st.just(300.0))
def test_case4_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    
    '''If the temperature of the left thermostat is equal to the temperature of the well and of the bar, 
    and it is higher than the temperature of the right thermostat, the temperature profile will be constant
    from the left thermostat to the well, and than decreasing till the right thermostat.'''
    
    bar_well = fp.bar_builder_well(temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar_well = fp.DFF_well(lenght,time,np.copy(bar_well),well_position)
    well_position=int(well_position/lenght*bar_well.shape[0])
    assert all(final_bar_well[i,-1] == final_bar_well[i+1,-1] for i in range(bar_well.shape[0]-(bar_well.shape[0]-well_position)))
    assert all(final_bar_well[well_position+i,-1] >= final_bar_well[well_position+i+1,-1] for i in range(bar_well.shape[0]-well_position-1))





