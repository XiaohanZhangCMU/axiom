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
    #if epoch > 60:
    #    return 1e-4
    #elif epoch > 100:
    #    return 1e-5
    return learning_rate


def solve_pde(model, **args):

    gstep = 0 # For restore model ckpts
    convg_tol = 1e-5 # convergence tol for float32
    epochs = args['epochs']
    LOG_DIR = args['LOG_DIR']

    GD, qdmr = dataset(reuse = False, **args)

    Gs = [ ]
    for Gid in sorted(GD.keys()):

        sess = tf.Session()
        sess.run(tf.global_variables_initializer())
        sess.run(tf.local_variables_initializer())

        if args['restore']:
            print('Restore saved model for reuse')
            saver = tf.train.Saver()
            saver.restore(sess, save_path=os.path.join(LOG_DIR, "Gid_"+str(Gid)+"_model.ckpt-0"))

        db = GD[Gid]
        inds = np.arange(db.shape[0])

        converged = False
        for epoch in range(epochs):
            np.random.shuffle(inds)

            for index, ind in enumerate(inds):
                result = sess.run([model.op, model.loss_sum, model.G], feed_dict={model.x_ph:db[ind], model.learning_rate:lr(epoch)})

            if epoch % 1 == 0:
                print("Epoch = {:5d}; LR = {:5.10E}; Residuals={: 5.10E};".format(epoch, lr(epoch), result[1]))

            if np.abs(result[1]) < convg_tol:
                cprint("Converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
                print("Residual={: 5.10E}.".format(result[1]))
                converged = True
                break
        if not converged:
            cprint("Failed to converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')

        saver = tf.train.Saver()
        saver.save(sess, os.path.join(LOG_DIR, "Gid_"+str(Gid)+"_model.ckpt"), gstep)



