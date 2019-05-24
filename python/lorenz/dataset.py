import os
import numpy as np

def dataset_sphere(reuse=False, **args):
    """ Generate dataset for space sampling
    """
    R = args['R']
    R0 = args['R0']
    center = args['center']
    n_samples = args['n_samples']
    batch_size = args['batch_size']

    db_file = args['db_save_path'] + 'open_cover_sphere_R_{}_Cx_{}_Cz_{}_batchsize_{}_nsamples_{}.npy'.format(R, center[0], center[1], batch_size, n_samples)

    print('R= {0}, Cx = {1}, Cz = {2}, batchsize = {3}'.format(R, center[0], center[1], batch_size))

    # Check if database is reusable
    if reuse and os.path.exists(db_file):
        return np.load(db_file)

    def sample_domain(batch_size=32):
        r  = np.sqrt(np.random.uniform(0, R, batch_size))
        angle = np.pi * np.random.uniform(0, 2, batch_size)
        xs = (r * np.cos(angle)+center[0]).reshape(-1,1)
        zs = (r * np.sin(angle)+center[1]).reshape(-1,1)
        return np.concatenate([xs,zs], axis=1)

    db = np.array([sample_domain(batch_size) for _ in range(n_samples)])
    np.save(db_file, db)
    return db

def dataset(reuse=False, **args):
    """ Generate dataset for space sampling
    """
    xmin, xmax = args['xmin'], args['xmax']
    zmin, zmax = args['zmin'], args['zmax']
    sub_nx, sub_nz = args['sub_nx'], args['sub_nz']
    patch_x, patch_z = args['patch_x'], args['patch_z']
    n_samples = args['n_samples']
    batch_size = args['batch_size']

    db_file = args['db_save_path'] + 'open_cover_xmin_{:5.2f}_xmax_{:5.2f}_zmin_{:.2f}_zmax_{:.2f}_subnx_{}_subnz_{}_patch_x_{:.2f}_patch_z_{:.2f}_batchsize_{}_nsamples_{}.npy'.format(xmin, xmax, zmin, zmax, sub_nx, sub_nz, patch_x, patch_z, batch_size, n_samples)

    print('xmin = {0}, xmax = {1}, zmin = {2}, zmax = {3}'.format(xmin, xmax, zmin, zmax))

    # Check if database is reusable
    if reuse and os.path.exists(db_file):
        return np.load(db_file)

    def sample_domain(batch_size=32):
        xs = np.random.uniform(xmin, xmax, batch_size).reshape(-1,1)
        zs = np.random.uniform(zmin, zmax, batch_size).reshape(-1,1)
        return np.concatenate([xs,zs], axis=1)

    db = np.array([sample_domain(batch_size) for _ in range(n_samples)])
    np.save(db_file, db)
    return db


