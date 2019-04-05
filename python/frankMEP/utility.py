""" A set of utility functions to process set operations of lists of nucleus indices
"""
import numpy as np
import pickle
#import matplotlib.pyplot as plt

def removeNucleusAtoms(cfg, nucleus, NP):
    idx = np.setdiff1d(np.xarrange(NP), nucleus)
    cfg = cfg[idx,:]
    return cfg

# Return a list of atom adjacent to, but not in, nucleus
def find_nbr_atoms(nbrlist, nucleus, totIdx):
    B = nbrlist[nucleus,:]
    A = np.setdiff1d(np.extract(B>=0,B), nucleus)
    return np.intersect1d(totIdx, A)

def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)

def merge_obj(obj_1, obj_2):
    obj = obj_1.copy()
    obj.update(obj_2)
    return obj

def nucleus2bits(nucleus, totIdx):
    bits = np.in1d(totIdx, nucleus)
    bits = bits.astype(int)
    return bits

def bits2nucleus(bits, totIdx):
    nucleus = totIdx[np.where(bits==1)]
    nucleus = nucleus.astype(int)
    return nucleus

def bits2str(bits):
    mybits=""
    for k in bits:
      mybits+=str(k)
    return mybits

def str2bits(bitstr):
    return np.array(map(int, bitstr))

def getnbrlist(nucleus, nbrlist):
    tmp = nbrlist[nucleus,:]
    return  np.extract(tmp>=0,tmp)

def nbitsdiff(bit1, bit2):
    #  print('diff by {0}'.format(len((np.where((bit1^bit2)==1))[0])))
    return len((np.where((bit1^bit2)==1))[0])

def nbits1(bit):
    return len(np.where(bit==1)[0])

def merge_database_files():
    obj1_file = 'DBFILES/db_0.05_1'
    obj2_file = 'DBFILES/db_0.05_2'
    obj3_file = 'DBFILES/db_0.05_3'
    obj4_file = 'DBFILES/db_0.05_4'
    obj5_file = 'DBFILES/db_0.05_5'

    obj12 = merge_obj(load_obj(obj1_file), load_obj(obj2_file));
    print("obj12 has {0} keys".format(len(obj12.keys())))
    obj34 = merge_obj(load_obj(obj3_file), load_obj(obj4_file));
    print("obj34 has {0} keys".format(len(obj34.keys())))
    obj1234 = merge_obj(obj12, obj34);
    obj = merge_obj(obj1234, load_obj(obj5_file));

    print("merged obj has {0} keys".format(len(obj.keys())))
    save_obj(obj, 'DBFILES/db_0.05');

def dbfile_sanity_check():
    db = load_obj('DBFILES/db_0.05_clearance');
    import numpy as np
    import matplotlib.pyplot as plt
    dbarr = []
    error_cnt = 0
    print('db has {0} keys.'.format(len(db.keys())))
    for key, val in db.items():
        dbarr.append([key.count('1'), val])
        if val > -165000:
            error_cnt += 1
            print('val = {0}, error count = {1}'.format(val, error_cnt))

    dbarr = np.array(dbarr)

    #plt.scatter(dbarr[:,0], dbarr[:,1])
    #plt.show()

def dbfile_rm_outliers():
    db = load_obj('DBFILES/db_0.05');
    print('db has {0} keys before clearance.'.format(len(db.keys())))
    for key in list(db.keys()):
        if db[key] > -165000:
            del db[key]

    print('db has {0} keys after clearance.'.format(len(db.keys())))
    save_obj(db, 'DBFILES/db_0.05_clearance');

if __name__ == "__main__":
    #merge_database_files()
    dbfile_rm_outliers()
    dbfile_sanity_check()

