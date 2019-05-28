import os
import pickle
import numpy as np

class QuadMapper:
    def __init__(self, nEx, nEy, domain):
        assert(len(domain)==4), "domain := [xmin, xmax, zmin, zmax]"
        self.quadMesher(nEx, nEy, domain)
        self.GD = {}
        self.build_GE() # Gid --> [ Eids ]
        self.build_EG() # Eid --> [ Gids ]
        # GD := Gid --> [ db file ]
        self.GD = {}
        for ky, val in self.GE.items():
            self.GD[ky] = 'db_'+str(ky) +'.npy'

    def build_GE(self):
        self.GE = {}
        id = 0
        for j in range(self.nEy-1):
            for i in range(self.nEx-1):
                eid1 = i + j * self.nEx
                eid2 = i + (j+1) * self.nEx
                self.GE[id] = [eid1, eid1+1, eid2, eid2+1]
                id += 1

    def build_EG(self):
        self.EG = {}
        for ky, val in self.GE.items():
            for eid in val:
                if eid not in self.EG:
                    self.EG[eid] = [ky]
                else:
                    self.EG[eid].append(ky)

    def locate_Gs(self, x):
        assert(len(x)==2)
        assert(x[0]<=self.domain[1] and x[0]>=self.domain[0] and x[1]<=self.domain[3] and x[1]>=self.domain[2])
        return self.EG[(x[0]-self.pad_domain[0])//self.dLx + (x[1]-self.pad_domain[2])//self.dLy * self.nEx]

    def quadMesher(self, nEx, nEy, domain):
        assert(nEx >= 1 and nEy >= 1)

        Lx0, Ly0 = domain[1] - domain[0] , domain[3] - domain[2]
        dLx, dLy = Lx0/nEx, Ly0/nEy

        # Padding with additional layer around domain
        pad_domain=[domain[0]-dLx, domain[1]+dLx, domain[2]-dLy, domain[3]+dLy]
        pad_Lx, pad_Ly=pad_domain[1] - pad_domain[0], pad_domain[3] - pad_domain[2]
        nEx, nEy = nEx + 2, nEy +2
        nnelem, nNodes, nElements = 4, (nEx+1)*(nEy+1), nEx*nEy

        rn0 = np.zeros((nNodes,2)) # Nodes ...
        for j in range(nEy+1):
            for i in range(nEx+1):
                rn0[i+j*(nEx+1),:] = [ (i/nEx-0.5)*pad_Lx + pad_domain[0],
                                       (j/nEy-0.5)*pad_Ly + pad_domain[2] ]

        elements = [[ 0 for _ in range(nnelem)] for _ in range(nElements)] # Elements ...
        for j in range(nEy):
            for i in range(nEx):
                elements[i+j*nEx] =[ a+i+j*(nEx+1) for a in [0, 1, nEx+2, nEx+1]];

        self.nEx, self.nEy = nEx, nEy
        self.rn0 = rn0
        self.elements = elements
        self.dLx, self.dLy = dLx, dLy
        self.domain = domain
        self.pad_domain = pad_domain

    def plot(self, X, x, conn, plotseq = [0,1,2,3,0]):
        import matplotlib.pyplot as plt
        dim, nnelem = X.shape[1], len(conn[0])
        for e in conn: plt.plot(X[e][plotseq,0],X[e][plotseq,1],'k:',x[e][plotseq,0],x[e][plotseq,1],'b-');
        plt.show()


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
    domain = [ args['xmin'], args['xmax'], args['zmin'], args['zmax'] ]
    nEx, nEy = args['Nx'], args['Ny']
    batch_size = args['batch_size']
    n_samples = args['n_samples']

    db_file = args['db_save_path'] + 'nEx_{}_nEy_{}_bsz_{}_nsmpl_{}.npy'.format(nEx, nEy, batch_size, n_samples)

    # Check if database is reusable
    if reuse and os.path.exists(db_file):
        return np.load(db_file)

    def sample_domain(xmin, xmax, zmin, zmax, batch_size=32):
        xs = np.random.uniform(xmin, xmax, batch_size).reshape(-1,1)
        zs = np.random.uniform(zmin, zmax, batch_size).reshape(-1,1)
        return np.concatenate([xs,zs], axis=1)

    qdmr = QuadMapper(nEx, nEy, domain)
    GE, elements, rn0 = qdmr.GE, qdmr.elements, qdmr.rn0

    GD = {}
    for ky, eids in GE.items():
        ids = []
        for e in eids:
            ids += elements[e]
        xmin, zmin = np.min(rn0[ids],axis=0)
        xmax, zmax = np.max(rn0[ids],axis=0)
        GD[ky] = np.array([sample_domain(xmin, xmax, zmin, zmax) for _ in range(n_samples)])
        #print('GD[{}].shape = {}'.format(ky, GD[ky].shape))

    pickle.dump(GD, open(db_file, 'wb'))
    return GD, qdmr

if __name__ == '__main__':
    qdmr = QuadMapper(2,2, [-20, 20, 0,40])
    print("QuadMapper info:")
    print('nEx = {}, nEy = {}, domain = {}'.format(qdmr.nEx, qdmr.nEy, qdmr.domain))
    print('elements = {}'.format(qdmr.elements))

    print('GD = \n')
    print(qdmr.GD)

    print('EG = \n')
    print(qdmr.EG)

    print('GE = \n')
    print(qdmr.GE)

    print('locate G = \n')
    print(qdmr.locate_Gs(np.array([-20,40])))

