from termcolor import colored, cprint
import numpy as np
import tensorflow as tf
from pde_model import Lorenz_Model
from dataset import dataset
from utils import *

def run_fine_lorenz(icpt, N = 10, **configs):
    import sys
    sys.path.append('../../lib/') #axiom/lib contains mdsw
    import lorenz
    lrz = lorenz.Sim()
    icpt = [-10,5,23]
    lrz.run_fine_ODE(N, configs['sigma'], configs['beta'], configs['gamma'], 0, icpt[0], icpt[1], icpt[2])
    traj = lrz.trajectory()
    return traj

def run_coarse_lorenz(icpt = [-10, 23], N=10, dt = 0.001, **configs):
    model = Lorenz_Model(**configs)
    sess = tf.Session()
    saver = tf.train.Saver()
    coarse_traj = np.zeros((N,3))

    # Scale x, y, z with gamma
    xi, zi = icpt[0]/configs['gamma'], icpt[1]/configs['gamma']

    GD, qdmr = dataset(reuse = False, **configs)

    sigma, beta, gamma = configs['sigma'], configs['beta'], configs['gamma']

    for ti in range(N):
        coarse_traj[ti] = [xi, 0, zi]

        Gids = qdmr.locate_Gs(np.array([xi,zi]))

        yi = 0
        for Gid in Gids:
            saver.restore(sess, save_path=os.path.join(configs['LOG_DIR'], "Gid_"+str(Gid)+"_model.ckpt-0"))
            test_input = np.array([xi,zi]).reshape(1,-1)
            loss_sum_val = sess.run([model.loss_sum], feed_dict={model.x_ph:test_input})
            print('Gid = {}, test output of loss_sum = {}'.format(Gid, loss_sum_val))
            y0 = sess.run([model.G], feed_dict={model.x_ph:test_input})
            yi += y0[0][0][0]
        yi /= len(Gids)

        xdot = sigma*(yi-xi)
        zdot = gamma*xi*yi - beta*zi
        xi = xi + xdot * dt
        zi = zi + zdot * dt

        print('ti = {}, xi = {}, zi = {}'.format(ti, xi, zi))

    return coarse_traj
