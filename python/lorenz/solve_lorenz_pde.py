"""
    Use tensorflow to solve PDE of lorenz-G equation
    1) Formulate least square loss
    2) Sample PDE domain patchwise

    Questions to resolve:
    1) Normalize input patch to have values [-1,1] for each channel?
    2) Add laplacian regularizer to PDE.
    3) Tune Adam, try other SGD optimizers, tune learning rate
    4) Try deeper neural network
    5) Initialize weights smarter ways
    6) Parallel speeding up GPU or MPI
    7) Prepare dataset beforehand, saving sampling time. Visualize each datapoint
    8) Add adaptive sampling scheme to sample more patches from high loss region

    Found papers:
    1) https://arxiv.org/pdf/1708.07469.pdf essentially following the same steps
"""
import sys
from termcolor import colored, cprint
import numpy as np
import tensorflow as tf
#from model import model
from dataset import dataset
from utils import *

def lr(epoch):
    learning_rate = 1e-3
    if epoch > 80:
        learning_rate *= 0.5e-3
    elif epoch > 60:
        learning_rate *= 1e-3
    elif epoch > 40:
        learning_rate *= 1e-2
    elif epoch > 20:
        learning_rate *= 1e-1
    return learning_rate


def solve_lorenz_pde(**args):

    tol = 1e-9 # convergence tol
    batch_size = args['batch_size'] # 1<<5
    epochs = args['epochs'] # 1<<20
    sigma = args['sigma']
    beta = args['beta']
    gamma = args['gamma']

    x_ph = tf.placeholder(dtype=np.float64, shape=(None, args['sub_nx'], args['sub_nz'], 2), name='input')
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
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate, beta1=0.9,beta2=0.999,epsilon=1e-08)
    grads_and_vars = optimizer.compute_gradients(loss_sum)
    op = optimizer.apply_gradients(grads_and_vars,global_step=tf.Variable(0, trainable=False))

    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())

    refine_factor = 0.8
    for refine_level in range(1):

        args['patch_x'] = args['patch_x']*refine_factor if args['patch_x']*refine_factor > args['min_patch_x'] else args['min_patch_x']
        args['patch_z'] = args['patch_z']*refine_factor if args['patch_z']*refine_factor > args['min_patch_z'] else args['min_patch_z']

        for epoch in range(epochs): # Minimize loss_sum til |loss-prev_loss|<tol

            db = dataset(reuse = False, **args)
            inds = np.arange(db.shape[0])
            np.random.shuffle(inds)

            for index, ind in enumerate(inds):
                result = sess.run([op, loss_sum, G], feed_dict={x_ph:db[ind], learning_rate:lr(epoch)})
                if index % 100 == 0:
                    print("Index = {:5d}; Residuals={: 5.10E};".format(index, result[1]))

            if epoch % 1 == 0:
                print("Epoch = {:5d}; LR = {:5.10E}; Residuals={: 5.10E};".format(epoch, lr(epoch), result[1]))
                # Save model weights
                # ...
            if np.abs(result[1]) < 1e-9:
                cprint("Converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
                print("Residual={: 5.10E}.".format(result[1], nrm, err))
                break

        cprint("Failed to converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')


#    # Test
#    xi = 10
#    zi = 23
#    N = 10
#    dt = 0.001
#
#    t =0
#    coarse_traj = np.zeros((N,3))
#    coarse_traj[t,0] = xi
#    coarse_traj[t,1] = 0
#    coarse_traj[t,2] = zi
#
#    for t in range(N):
#        test_input = np.array(np.stack(np.meshgrid(np.linspace(xi-.5*patch_x,xi+.5*patch_x,sub_nx),
#        np.linspace(zi-.5*patch_z, zi+.5*patch_z, sub_nz)), axis=-1))
#        new_test_input = test_input.reshape((1, *(test_input.shape)))
#        loss_sum_val = sess.run([loss_sum], feed_dict={x_ph:new_test_input})
#        print('test output of loss_sum = {0}'.format(loss_sum_val[0]))
#        y0 = sess.run([G], feed_dict={x_ph:new_test_input})
#        yi = y0[0][0,7,7,0]
#        xdot = sigma*(yi-xi)
#        zdot = xi*yi - beta*zi
#        xi = xi + xdot * dt
#        zi = zi + zdot * dt
#
#        coarse_traj[t,0] = xi
#        coarse_traj[t,1] = 0
#        coarse_traj[t,2] = zi

    step = 0
    saver = tf.train.Saver()
    LOG_DIR = "/Users/xiaohan.zhang/Planet/Codes/axiom/python/lorenz/"
    saver.save(sess, os.path.join(LOG_DIR, "model.ckpt"), step)

    # Load from saved model
    # saver = tf.train.Saver()
    # saver.restore(sess, save_path=os.path.join(LOG_DIR, "model.ckpt-10"))

if __name__ == '__main__':

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

    result= solve_lorenz_pde(n_samples=n_samples, batch_size=batch_size,
            xmin=xmin, xmax=xmax, zmin=zmin, zmax=zmax,
            sub_nx=sub_nx, sub_nz=sub_nz,
            patch_x=patch_x, patch_z=patch_z,
            min_patch_x = min_patch, min_patch_z = min_patch,
            db_save_path=db_save_path,
            sigma=sigma, beta=beta, gamma=gamma, epochs=epochs)






