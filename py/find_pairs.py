from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
import numpy as np
from sklearn.neighbors import NearestNeighbors
import h5py

def load_abacus_summit(sim_path):
    print()
    print(sim_path)
    print("Started reading the data")
    cat = CompaSOHaloCatalog(sim_path, fields=["vcirc_max_L2com", "x_L2com", "v_L2com"])
    print("Finished reading the data")
    
    n_halos = len(cat.halos)
    print('Loaded {:.1f}M halos in total'.format(n_halos/1E6))
    
    # filtering only massive halos
    cat.halos['ID'] = np.arange(len(cat.halos))
    ii = cat.halos['vcirc_max_L2com']>220.0
    halo_data = cat.halos[ii]
    n_halos = len(halo_data)
    print('Kept {:.1f}M halos in total'.format(n_halos/1E6))

    # sorting
    print("Started sorting the data")
    halo_data = np.sort(halo_data, order='vcirc_max_L2com')
    print("Finished sorting the data")
    
    S_pos = halo_data['x_L2com']
    S_vel = halo_data['v_L2com']
    S_vmax = halo_data['vcirc_max_L2com']
    S_mass = halo_data['N']
    S_parent_fof = halo_data['ID']
    
    return S_pos, S_vel, S_vmax, S_mass, S_parent_fof

def find_pairs(S_pos, S_vel, S_vmax, S_mass, S_parent_fof, pair_filename):
    n_S = len(S_pos)
    print("Number of halos to use:", n_S)
    
    print("Started Neighbor computation")
    nbrs_S = NearestNeighbors(n_neighbors=20, algorithm='ball_tree').fit(S_pos)
    dist_S, ind_S = nbrs_S.kneighbors(S_pos)
    print("Finished Neighbor computation")
    
    print("Started isolation computation")
    neighbor_index = ind_S[:,1]
    neighbor_list = ind_S[:,2:]

    n_pairs = 0

    halo_A_id = np.empty((0), dtype=int)
    halo_B_id = np.empty((0), dtype=int)

    for i in range(n_S):
        l = neighbor_index[neighbor_index[i]]% n_S
        j = neighbor_index[i] % n_S
    
        other_j = neighbor_list[i,:] % n_S
        other_l = neighbor_list[neighbor_index[i],:] % n_S
    
        if((i==l) & (not (j in halo_A_id)) & (not (i in halo_B_id))): # first check to find mutual neighbors
            vmax_i = S_vmax[i]
            vmax_j = S_vmax[j]
            vmax_limit = min([vmax_i, vmax_j])
                
            pair_d = dist_S[i,1] # This is the current pair distance
            dist_limit = pair_d * 3.0 # exclusion radius for massive structures
            
            massive_close_to_i = any((dist_S[i,2:]<dist_limit) & (S_vmax[other_j] >= vmax_limit))
            massive_close_to_j = any((dist_S[j,2:]<dist_limit) & (S_vmax[other_l] >= vmax_limit))
            if((not massive_close_to_i) & (not massive_close_to_j)): # check on massive structures inside exclusion radius
                n_pairs = n_pairs+ 1
                halo_A_id = np.append(halo_A_id, int(i))
                halo_B_id = np.append(halo_B_id, int(j))
    print("Finished isolation computation")
    print("Pairs found:", n_pairs)
    
    
    print("Started writing data to ", pair_filename)
    h5f = h5py.File(pair_filename, 'w')
    h5f.create_dataset('pos_A', data=S_pos[halo_A_id,:])
    h5f.create_dataset('pos_B', data=S_pos[halo_B_id,:])
    h5f.create_dataset('N_A', data=S_mass[halo_A_id])
    h5f.create_dataset('N_B', data=S_mass[halo_B_id])
    #h5f.create_dataset('pos_G', data=S_pos)
    h5f.create_dataset('vel_A', data=S_vel[halo_A_id,:])
    h5f.create_dataset('vel_B', data=S_vel[halo_B_id,:])
    #h5f.create_dataset('vel_G', data=S_vel)
    h5f.create_dataset('vmax_A', data=S_vmax[halo_A_id])
    h5f.create_dataset('vmax_B', data=S_vmax[halo_B_id])
    #h5f.create_dataset('vmax_G', data=S_vmax)
    return 

# 2 Gpc/h, N=6912^3 particles, base resulution. Fixed cosmology, changing phases
#for phase in range(5):
#    sim_path = "/global/cfs/cdirs/desi/cosmosim/Abacus/AbacusSummit_base_c000_ph{:03d}/halos/z0.100".format(phase)
#    S_pos, S_vel, S_vmax, S_mass, S_parent_fof = load_abacus_summit(sim_path)
#    pair_filename = "../data/pairs_AbacusSummit_base_c000_ph{:03d}_z0.100.hdf5".format(phase)
#    find_pairs(S_pos, S_vel, S_vmax, S_mass, S_parent_fof, pair_filename)
    
# 2 Gpc/h, N=6912^3 particles, base resulution. Fixed phase, different cosmologies
# ver tabla 3 https://arxiv.org/pdf/2110.11398.pdf
# Illustris es la cosmologia 017
#for cosmo in range(14,19):
#    sim_path = "/global/cfs/cdirs/desi/cosmosim/Abacus/AbacusSummit_base_c{:03d}_ph000/halos/z0.100".format(cosmo)
#    S_pos, S_vel, S_vmax, S_mass, S_parent_fof = load_abacus_summit(sim_path)
#    pair_filename = "../data/pairs_AbacusSummit_base_c{:03d}_ph000_z0.100.hdf5".format(cosmo)
#    find_pairs(S_pos, S_vel, S_vmax, S_mass, S_parent_fof, pair_filename)

# 1 Gpc/h, base mass resolution. 
#sim_path = "/global/cfs/cdirs/desi/cosmosim/Abacus/AbacusSummit_highbase_c000_ph100/halos/z0.100"
#S_pos, S_vel, S_vmax, S_mass, S_parent_fof = load_abacus_summit(sim_path)
#pair_filename = "../data/pairs_AbacusSummit_highbase_c000_ph100_z0.100.hdf5"
#find_pairs(S_pos, S_vel, S_vmax, S_mass, S_parent_fof, pair_filename)


# 1 Gpc/h, N=6300^3 particles, high resolution.
#sim_path = "/global/cfs/cdirs/desi/cosmosim/Abacus/AbacusSummit_high_c000_ph100/halos/z0.100"
#S_pos, S_vel, S_vmax, S_mass, S_parent_fof = load_abacus_summit(sim_path)
#pair_filename = "../data/pairs_AbacusSummit_high_c000_ph100_z0.100.hdf5"
#find_pairs(S_pos, S_vel, S_vmax, S_mass, S_parent_fof, pair_filename)

#sim_path = "/global/cfs/cdirs/desi/cosmosim/Abacus/AbacusSummit_high_c000_ph100/halos/z0.100"
#S_pos, S_vel, S_vmax, S_mass, S_parent_fof = load_abacus_summit(sim_path)
#pair_filename = "../data/pairs_AbacusSummit_high_c000_ph100_z0.100.hdf5"
#find_pairs(S_pos, S_vel, S_vmax, S_mass, S_parent_fof, pair_filename)


    


