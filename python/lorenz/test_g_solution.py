import sys
from termcolor import colored, cprint
import numpy as np
import tensorflow as tf
#from model import model
from dataset import dataset
from utils import *

n_samples = 1<<10
batch_size = 1<<5
sub_nx, sub_nz = 16, 16
patch_x, patch_z = 2.5, 2.5
min_patch = 0.025
xmin, xmax = -20, 20
zmin, zmax = 0, 40
db_save_path = './'
epsilon = 0.1
sigma = 10
beta = 8./3
gamma = 25
epochs=2

x_ph = tf.placeholder(dtype=np.float64, shape=(None, sub_nx, sub_nz, 2), name='input')
learning_rate = tf.placeholder(tf.float32, shape=[], name='learning_rate')

def mlp(x, hidden_sizes=(32,1), activation=tf.tanh, output_activation=None):
    for h in hidden_sizes[:-1]:
        x = tf.layers.dense(x, units=h, activation=activation)
    return tf.layers.dense(x, units=hidden_sizes[-1], activation=output_activation)

def loss_function(inputs):
    xadv = tf.identity(inputs) # Such that dG/d_xadv = dG/d_x
    G = mlp(xadv, hidden_sizes=(16,1))
    dG = tf.gradients(G, xadv)[0]
    # d2G,_ = tf.hessians(G, xadv) # Regularizer
    print('xadv, G and dG shape')
    print(xadv.shape)
    print(G.shape)
    print(dG.shape)
    print((G-xadv[:,:,:,0:1]).shape)

    loss_sum = tf.reduce_mean(((sigma*((G-xadv[:,:,:,0:1])*dG[:,:,:,0:1])) +
        (((xadv[:,:,:,0:1]*G)-beta*xadv[:,:,:,1:2])*dG[:,:,:,1:2]) +
                               G +
                               (xadv[:,:,:,0:1]*(xadv[:,:,:,1:2]-gamma))
                                # + tf.math.scalar_mul(epsilon,(d2Gx+d2Gz))  # Regularizer
                               )**2)
    print(loss_sum.shape)
    print('done')
    return loss_sum, G

    """ Optimization setup
    """

loss_sum,G = loss_function(x_ph)
sess = tf.Session()
step = 0
LOG_DIR = "/Users/xiaohan.zhang/Planet/Codes/axiom/python/lorenz/"
saver = tf.train.Saver()
saver.restore(sess, save_path=os.path.join(LOG_DIR, "model.ckpt-0"))

# Test
xi = 10
zi = 23
N = 10
dt = 0.001

t =0
coarse_traj = np.zeros((N,3))
coarse_traj[t,0] = xi
coarse_traj[t,1] = 0
coarse_traj[t,2] = zi

for t in range(N):
    test_input = np.array(np.stack(np.meshgrid(np.linspace(xi-.5*patch_x,xi+.5*patch_x,sub_nx),
    np.linspace(zi-.5*patch_z, zi+.5*patch_z, sub_nz)), axis=-1))
    new_test_input = test_input.reshape((1, *(test_input.shape)))
    loss_sum_val = sess.run([loss_sum], feed_dict={x_ph:new_test_input})
    print('test output of loss_sum = {0}'.format(loss_sum_val[0]))
    y0 = sess.run([G], feed_dict={x_ph:new_test_input})
    yi = y0[0][0,7,7,0]
    xdot = sigma*(yi-xi)
    zdot = xi*yi - beta*zi
    xi = xi + xdot * dt
    zi = zi + zdot * dt

    coarse_traj[t,0] = xi
    coarse_traj[t,1] = 0
    coarse_traj[t,2] = zi

