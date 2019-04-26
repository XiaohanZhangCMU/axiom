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
    saver.restore(sess, save_path=os.path.join(configs['LOG_DIR'], "model.ckpt-0"))
    coarse_traj = np.zeros((N,3))
    xi, zi = icpt[0], icpt[1]

    min_patch_x = configs['min_patch_x']
    min_patch_z = configs['min_patch_z']
    sub_nx = configs['sub_nx']
    sub_nz = configs['sub_nz']
    sigma, beta, gamma = configs['sigma'], configs['beta'], configs['gamma']

    for ti in range(N):
        coarse_traj[ti] = [xi, 0, zi]
        test_input = np.array(np.stack(np.meshgrid(np.linspace(xi-.5*min_patch_x,xi+.5*min_patch_x,sub_nx),np.linspace(zi-.5*min_patch_z, zi+.5*min_patch_z, sub_nz)), axis=-1))
        new_test_input = test_input.reshape((1, *(test_input.shape)))
        loss_sum_val = sess.run([model.loss_sum], feed_dict={model.x_ph:new_test_input})
        print('test output of loss_sum = {0}'.format(loss_sum_val[0]))
        y0 = sess.run([model.G], feed_dict={model.x_ph:new_test_input})
        yi = y0[0][0,7,7,0]
        xdot = sigma*(yi-xi)
        zdot = xi*yi - beta*zi
        xi = xi + xdot * dt
        zi = zi + zdot * dt

    return coarse_traj
