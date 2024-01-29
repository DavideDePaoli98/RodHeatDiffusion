from DiffEqLibrary import DiffEqLibrary_ as fp
import numpy as np

# Test to control if the return of the bar builder have the desired dimension and have the fixed temperatures in the extremetes 
# (and in the well, if it exists) during the entire simulation time
def test_bar_builder(temperature_left,temperature_right,temperature_bar):
    bar=fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    assert all(i == temperature_right for i in bar[0,:])
    assert all(i == temperature_left for i in bar[-1,:])
    assert all(i == temperature_bar for i in bar[1:-1,0]) 
def test_bar_builder_well(temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght):
    bar=fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    well_position=int(well_position/lenght*bar.shape[0])
    assert all(i == temperature_right for i in bar[0,:])
    assert all(i == temperature_left for i in bar[-1,:])
    assert all(i == temperature_well for i in bar[well_position,:])
    assert all(i == temperature_bar for i in bar[1:well_position,0])
    assert all(i == temperature_bar for i in bar[well_position+1:-1,0])

# Test to control if the return of the finite difference method functions have the desired dimension and have the fixed temperatures in the extremetes 
# (and in the well, if it exists) during the entire simulation time. It also controls if the return is different from the initial bar configuration.
def test_DFF(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.DFF(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_bar for i in final_bar[1:-1,0])
    assert not np.array_equal(initial_bar,final_bar)
def test_DFF_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.DFF_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])
    assert all(i == temperature_bar for i in final_bar[1:well_position,0])
    assert all(i == temperature_bar for i in final_bar[well_position+1:-1,0])
    assert not np.array_equal(initial_bar,final_bar)
def test_C_N(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.C_N(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert not np.array_equal(initial_bar,final_bar)
def test_C_N_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.C_N_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])
    assert not np.array_equal(initial_bar,final_bar)
def test_R_K(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.R_K(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert not np.array_equal(initial_bar,final_bar)
def test_R_K_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.R_K_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])
    assert not np.array_equal(initial_bar,final_bar)
def test_FTCS(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left,temperature_right,temperature_bar)
    final_bar = fp.FTCS(lenght,time,np.copy(initial_bar))
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert not np.array_equal(initial_bar,final_bar)
def test_FTCS_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left,temperature_right,temperature_bar,well_position,temperature_well,lenght)
    final_bar = fp.FTCS_well(lenght,time,np.copy(initial_bar),well_position)
    well_position=int(well_position/lenght*initial_bar.shape[0])
    assert all(i == temperature_right for i in final_bar[0,:])
    assert all(i == temperature_left for i in final_bar[-1,:])
    assert all(i == temperature_well for i in final_bar[well_position,:])
    assert not np.array_equal(initial_bar,final_bar)

#Test to control if the different functions work well with general initial conditions
test_bar_builder(100,200,300)
test_bar_builder_well(100,200,300,0.5,250,1)
test_DFF(100,200,300,1,10)
test_DFF_well(100,200,300,1,10,0.5,250)
test_C_N(100,200,300,1,10)
test_C_N_well(100,200,300,1,10,0.5,250)
test_R_K(100,200,300,1,10)
test_R_K_well(100,200,300,1,10,0.5,250)
test_FTCS(100,200,300,1,10)
test_FTCS_well(100,200,300,1,10,0.5,250)