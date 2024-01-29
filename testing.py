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
    initial_bar = fp.bar_builder (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.DFF(lenght,time,initial_bar,linear_diffusion)
    assert final_bar.shape[0] == dim_X
    assert final_bar.shape[1] == dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1:
            assert final_bar[m][0] == temperature_bar
            assert final_bar[m][1] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]  == temperature_right
        assert final_bar[-1][n] == temperature_left
    assert not np.array_equal(initial_bar,final_bar)
def test_DFF_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,well_position_=well_position,temperature_well_=temperature_well,lenght_=lenght,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.DFF_well(lenght,time,initial_bar,well_position,linear_diffusion)
    well_position=int(well_position/lenght*dim_X)
    assert final_bar.shape[0]==dim_X
    assert final_bar.shape[1]==dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1 and m != well_position:
            assert final_bar[m][0] == temperature_bar
            assert final_bar[m][1] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]             == temperature_right
        assert final_bar[-1][n]            == temperature_left
        assert final_bar[well_position][n] == temperature_well
    assert not np.array_equal(initial_bar,final_bar)
def test_C_N(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.C_N(lenght,time,initial_bar,linear_diffusion)
    assert final_bar.shape[0] == dim_X
    assert final_bar.shape[1] == dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1:
            assert final_bar[m][0] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]  == temperature_right
        assert final_bar[-1][n] == temperature_left
    assert not np.array_equal(initial_bar,final_bar)
def test_C_N_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,well_position_=well_position,temperature_well_=temperature_well,lenght_=lenght,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.C_N_well(lenght,time,initial_bar,well_position,linear_diffusion)
    well_position=int(well_position/lenght*dim_X)
    assert final_bar.shape[0]==dim_X
    assert final_bar.shape[1]==dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1 and m != well_position:
            assert final_bar[m][0] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]             == temperature_right
        assert final_bar[-1][n]            == temperature_left
        assert final_bar[well_position][n] == temperature_well
    assert not np.array_equal(initial_bar,final_bar)
def test_R_K(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.R_K(lenght,time,initial_bar,linear_diffusion)
    assert final_bar.shape[0] == dim_X
    assert final_bar.shape[1] == dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1:
            assert final_bar[m][0] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]  == temperature_right
        assert final_bar[-1][n] == temperature_left
    assert not np.array_equal(initial_bar,final_bar)
def test_R_K_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,well_position_=well_position,temperature_well_=temperature_well,lenght_=lenght,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.R_K_well(lenght,time,initial_bar,well_position,linear_diffusion)
    well_position=int(well_position/lenght*dim_X)
    assert final_bar.shape[0]==dim_X
    assert final_bar.shape[1]==dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1 and m != well_position:
            assert final_bar[m][0] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]             == temperature_right
        assert final_bar[-1][n]            == temperature_left
        assert final_bar[well_position][n] == temperature_well
    assert not np.array_equal(initial_bar,final_bar)
def test_FTCS(temperature_left,temperature_right,temperature_bar,lenght,time):
    initial_bar = fp.bar_builder (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.FTCS(lenght,time,initial_bar,linear_diffusion)
    assert final_bar.shape[0] == dim_X
    assert final_bar.shape[1] == dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1:
            assert final_bar[m][0] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]  == temperature_right
        assert final_bar[-1][n] == temperature_left
    assert not np.array_equal(initial_bar,final_bar)
def test_FTCS_well(temperature_left,temperature_right,temperature_bar,lenght,time,well_position,temperature_well):
    initial_bar = fp.bar_builder_well (temperature_left_=temperature_left,temperature_right_=temperature_right,temperature_bar_=temperature_bar,well_position_=well_position,temperature_well_=temperature_well,lenght_=lenght,dim_X_=dim_X,dim_t_=dim_t)
    final_bar = fp.FTCS_well(lenght,time,initial_bar,well_position,linear_diffusion)
    well_position=int(well_position/lenght*dim_X)
    assert final_bar.shape[0]==dim_X
    assert final_bar.shape[1]==dim_t
    for m in range(dim_X):
        if m != 0 and m != dim_X-1 and m != well_position:
            assert final_bar[m][0] == temperature_bar
    for n in range(dim_t):
        assert final_bar[0][n]             == temperature_right
        assert final_bar[-1][n]            == temperature_left
        assert final_bar[well_position][n] == temperature_well
    assert not np.array_equal(initial_bar,final_bar)

#Test to control the different functions returns considering variation about the lenght of the bar and the time interval of the simulation
for t in range (50):
    time=t+1
    for x in range(5):
        lenght=x+1
        test_bar_builder(0,10,50,100,100)
        test_bar_builder_well(0,10,50,lenght/2,10,lenght,100,100)
        test_DFF(0,10,50,lenght,time,0.00005,100,100)
        test_DFF_well(0,10,50,lenght,time,0.00005,lenght/2,10,100,100)
        test_C_N(0,10,50,lenght,time,0.00005,100,100)
        test_C_N_well(0,10,50,lenght,time,0.00005,lenght/2,10,50,50)
        test_R_K(0,10,50,lenght,time,0.00005,100,100)
        test_R_K_well(0,10,50,lenght,time,0.00005,lenght/2,10,100,100)
        test_FTCS(0,10,50,lenght,time,0.00005,100,100)
        test_FTCS_well(0,10,50,lenght,time,0.00005,lenght/2,10,100,100)
#Test to control the different functions returns considering variation about temperature of one of the different thermostat
for temp in range(10):
    lenght= 2
    time= 50
    temper = 10**temp
    test_bar_builder(0,10,temper,100,100)
    test_bar_builder_well(0,10,temper,lenght/2,10,lenght,100,100)
    test_DFF(0,10,temper,lenght,time,0.00005,100,100)
    test_DFF_well(0,10,50,lenght,time,0.00005,lenght/2,10,100,100)
    test_C_N(0,10,temper,lenght,time,0.00005,100,100)
    test_C_N_well(0,temper,50,lenght,time,0.00005,lenght/2,10,50,50)
    test_R_K(0,10,temper,lenght,time,0.00005,100,100)
    test_R_K_well(0,10,temper,lenght,time,0.00005,lenght/2,10,100,100)
    test_FTCS(0,10,50,lenght,time,0.00005,100,100)
    test_FTCS_well(0,10,temper,lenght,time,0.00005,lenght/2,10,100,100)
