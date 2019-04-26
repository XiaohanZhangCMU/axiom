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
from pde_model import Lorenz_Model
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


def solve_pde(model, **args):

    tol = 1e-9 # convergence tol
    epochs = args['epochs']
    LOG_DIR = args['LOG_DIR']

    """ Optimization setup
    """

    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())
    refine_factor = 0.8

    for refine_level in range(1):

        db = dataset(reuse = False, **args)
        inds = np.arange(db.shape[0])

        for epoch in range(epochs): # Minimize loss_sum til |loss-prev_loss|<tol

            np.random.shuffle(inds)

            for index, ind in enumerate(inds):
                result = sess.run([model.op, model.loss_sum, model.G], feed_dict={model.x_ph:db[ind], model.learning_rate:lr(epoch)})
                if index == 0 or index % 100 == 0:
                    print("Index = {:5d}; Residuals={: 5.10E};".format(index, result[1]))

            if epoch % 1 == 0:
                print("Epoch = {:5d}; LR = {:5.10E}; Residuals={: 5.10E};".format(epoch, lr(epoch), result[1]))

            if np.abs(result[1]) < 1e-9:
                cprint("Converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
                print("Residual={: 5.10E}.".format(result[1], nrm, err))
                break

        cprint("Failed to converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')

        args['patch_x'] = args['patch_x']*refine_factor if args['patch_x']*refine_factor > args['min_patch_x'] else args['min_patch_x']
        args['patch_z'] = args['patch_z']*refine_factor if args['patch_z']*refine_factor > args['min_patch_z'] else args['min_patch_z']



    step = 0
    saver = tf.train.Saver()
    saver.save(sess, os.path.join(LOG_DIR, "model.ckpt"), step)


