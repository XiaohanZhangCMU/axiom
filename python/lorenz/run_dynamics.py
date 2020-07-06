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
    lrz.run_fine_ODE(N, configs['sigma'], configs['beta'], configs['gamma'], 0, icpt[0], icpt[1], icpt[2])
    traj = lrz.trajectory()
    return traj

def run_coarse_lorenz(icpt = [-10, 23], N=10, dt = 0.001, **configs):
    model = Lorenz_Model(**configs)
    sess = tf.Session()
    saver = tf.train.Saver()
    coarse_traj = np.zeros((N,3))

    # Scale x, y, z with gamma
    sigma, beta, gamma, Bbar = configs['sigma'], configs['beta'], configs['gamma'], configs['Bbar']
    xi, zi = icpt[0]/Bbar, icpt[1]/Bbar

    GD, qdmr = dataset(reuse = False, **configs)


    for ti in range(N):
        coarse_traj[ti] = [xi, 0, zi]

        Gids = qdmr.locate_Gs(np.array([xi,zi]))

        loss_sums, yis = [], []
        for Gid in Gids:
            saver.restore(sess, save_path=os.path.join(configs['LOG_DIR'], "Gid_"+str(Gid)+"_model.ckpt-0"))
            test_input = np.array([xi,zi]).reshape(1,-1)
            loss_sum_val, y0 = sess.run([model.loss_sum, model.G], feed_dict={model.x_ph:test_input})
            print('Gid = {}, test output of loss_sum = {}, y0 = {}'.format(Gid, loss_sum_val, y0))
            #y0 = sess.run([model.G], feed_dict={model.x_ph:test_input})
            #yi += y0[0][0][0]
            yis.append(y0[0][0])
            loss_sums.append(loss_sum_val)
        #yi /= len(Gids)
        yi = yis[np.argmin(loss_sums)]
        print('yi = {}'.format(yi))

        xdot = sigma*(yi-xi)
        zdot = Bbar*xi*yi - beta*zi
        xi = xi + xdot * dt
        zi = zi + zdot * dt

        print('ti = {}, xi = {}, zi = {}'.format(ti, xi, zi))

    return coarse_traj
