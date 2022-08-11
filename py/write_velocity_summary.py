import numpy as np
import matplotlib.pyplot as plt
import h5py
import os


def write_velocity_summary(filename_in, filename_out, illustris=False):
    data = {}
    f = h5py.File(filename_in, 'r')
    for k in f.keys():
        data[k] = f[k][...]
    f.close()
    print('finished reading {}'.format(filename_in))


    a = data['vel_A'].copy()
    b = data['vel_B'].copy()
    v_cm = data['vel_A'].copy()
    p_cm = data['pos_A'].copy()
    mass_tot = data['N_A'] + data['N_B']
    for i in range(3):
        a[:,i] = data['vel_A'][:,i] * data['N_A']/mass_tot
        b[:,i] = data['vel_B'][:,i] * data['N_B']/mass_tot
        v_cm[:,i] = a[:,i] + b[:,i]
        
        a[:,i] = data['pos_A'][:,i] * data['N_A']/mass_tot
        b[:,i] = data['pos_B'][:,i] * data['N_B']/mass_tot
        p_cm[:,i] = a[:,i] + b[:,i]
        
    v_cm_norm = np.sqrt(np.sum(v_cm**2, axis=1))
    data['vel_A_mag'] = np.sqrt(np.sum(data['vel_A']**2, axis=1))

    data['vel_A_mag'] = np.sqrt(np.sum(data['vel_A']**2, axis=1))
    data['vel_B_mag'] = np.sqrt(np.sum(data['vel_B']**2, axis=1))

    data['pos_AB'] = np.sqrt(np.sum( (data['pos_B'] - data['pos_A'])**2, axis=1))
    data['vel_AB'] = np.sqrt(np.sum( (data['vel_B'] - data['vel_A'])**2, axis=1))
    data['vel_AB_rad'] = np.sum((data['pos_B'] - data['pos_A'])*(data['vel_B'] - data['vel_A']), axis=1)/data['pos_AB']
    data['vel_AB_tan'] = np.sqrt((data['vel_AB']**2 - data['vel_AB_rad']**2))
    data['pos_CM'] = p_cm

    #now we compute the radial velocity including the hubble flow
    if illustris==False:
        data['vel_AB_rad'] = data['vel_AB_rad'] + (data['pos_AB'] * 100.0)
        ii = (data['vel_AB_rad']<0) & (data['pos_AB']<1.0) & ((data['vmax_A']<260.0) & (data['vmax_B']<260.0)) 
        ii = ii & (data['N_A']>0) & (data['N_B']>0)
    else:
        data['vel_AB_rad'] = data['vel_AB_rad'] + ((data['pos_AB']/1000)*100.0)
        ii = (data['vel_AB_rad']<0) & ((data['pos_AB']/1000.0)<1.0) & ((data['vmax_A']<260.0) & (data['vmax_B']<260.0)) 
        ii = ii & (data['N_A']>0) & (data['N_B']>0)

        

    cm_vel = v_cm_norm[ii]
    tan_vel = data['vel_AB_tan'][ii]
    rad_vel = data['vel_AB_rad'][ii]
    tot_mass = data['N_A'][ii] + data['N_B'][ii]
    pos_cm = data['pos_CM'][ii]
    
    results  = np.array([cm_vel, tan_vel, rad_vel, tot_mass, pos_cm[:,0], pos_cm[:,1], pos_cm[:,2]])

    print(np.shape(results))

    # write positions
    np.savetxt(filename_out, results.T, fmt='%f %f %f %f %f %f %f', 
              header='v_cm v_tan, v_rad m_tot x_cm y_cm z_cm')
    print(' wrote results data to {}'.format(filename_out))


for box in [50,100,300]:
    for res in [1,2,3]:
        filename_in = '../data/pairs_TNG{}-{}.hdf5'.format(box,res)
        filename_out = '../data/lg_pairs_TNG{}-{}.dat'.format(box,res)
        write_velocity_summary(filename_in, filename_out, illustris=True)
        
for phase in [0,1]:
    filename_in = '../data/pairs_AbacusSummit_base_c000_ph{:03d}_z0.100.hdf5'.format(phase)
    filename_out = '../data/lg_pairs_AbacusSummit_base_c000_ph{:03d}_z0.100.dat'.format(phase)
    write_velocity_summary(filename_in, filename_out, illustris=False)