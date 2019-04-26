""" Generate dataset for space sampling
"""
#def dataset(n_samples = 1<<15, batch_size=1<<5,
#            xmin=-1, xmax=1, zmin=-1, zmax=1,
#            sub_nx=16, sub_nz=16, patch_x=1., patch_z=1.,
#            db_save_path='./'):

import numpy as np

def dataset(reuse=False, **args):
    import os
    # Check if database is reusable

    xmin, xmax = args['xmin'], args['xmax']
    zmin, zmax = args['zmin'], args['zmax']
    sub_nx, sub_nz = args['sub_nx'], args['sub_nz']
    patch_x, patch_z = args['patch_x'], args['patch_z']
    n_samples = args['n_samples']
    batch_size = args['batch_size']

    db_file = args['db_save_path'] + 'open_cover_xmin_{:5.2f}_xmax_{:5.2f}_zmin_{:.2f}_zmax_{:.2f}_subnx_{}_subnz_{}_patch_x_{:.2f}_patch_z_{:.2f}_batchsize_{}_nsamples_{}.npy'.format(xmin, xmax, zmin, zmax, sub_nx, sub_nz, patch_x, patch_z, batch_size, n_samples)

    if reuse and os.path.exists(db_file):
        db = np.load(db_file)
        return db

    def sample_domain(batch_size=32):
        xs = np.random.uniform(xmin, xmax, batch_size)
        zs = np.random.uniform(zmin, zmax, batch_size)
        return np.array([np.stack(np.meshgrid(np.linspace(xi-.5*patch_x,xi+.5*patch_x,sub_nx),
            np.linspace(zi-.5*patch_z, zi+.5*patch_z, sub_nz)), axis=-1)
            for (xi,zi) in zip(xs,zs)])

    db = []
    print('batch_size = {0}'.format(batch_size))
    for ind in range(n_samples):
        X_sample = sample_domain(batch_size=batch_size)
        db.append(X_sample)
        if ind % 1000 == 0:
            print('dataset generating .... {0}/{1}'.format(ind, n_samples))

    db = np.array(db)
    print('Database of shape: {0} saved to {1}'.format(db.shape, db_file))
    np.save(db_file, db)

    return db


