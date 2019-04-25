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


#def solve_lorenz_pde(n_samples, batch_size, 
#                     xmin, xmax, zmin, zmax, 
#                     sub_nx, sub_nz, patch_x, patch_z, 
#                     db_save_path,
#                     epsilon, sigma, beta, gamma):

def solve_lorenz_pde(**args):

    lr = args['lr'] # 1e-3  # learning rate
    tol = 1e-9 # convergence tol
    batch_size = args['batch_size'] # 1<<5
    # n_batch = args['n_batch'] # 1<<6
    epochs = args['epochs'] # 1<<20
    sigma = args['sigma']
    beta = args['beta']
    gamma = args['gamma']

    x_ph = tf.placeholder(dtype=np.float64, shape=(None, args['sub_nx'], args['sub_nz'], 2), name='input')

#    def sample_domain(batch_size=32):
#        xmin, xmax = args['xmin'], args['xmax']
#        zmin, zmax = args['zmin'], args['zmax']
#        sub_nx, sub_nz = args['sub_nx'], args['sub_nz']
#        patch_x, patch_z = args['patch_x'], args['patch_z'] 
#
#        xs = np.random.uniform(xmin, xmax, batch_size)
#        zs = np.random.uniform(zmin, zmax, batch_size)
#        return np.array([np.stack(np.meshgrid(np.linspace(xi-.5*patch_x,xi+.5*patch_x,sub_nx),
#            np.linspace(zi-.5*patch_z, zi+.5*patch_z, sub_nz)), axis=-1)
#            for (xi,zi) in zip(xs,zs)])

    def mlp(x, hidden_sizes=(32,1), activation=tf.tanh, output_activation=None):
        for h in hidden_sizes[:-1]:
            x = tf.layers.dense(x, units=h, activation=activation)
        return tf.layers.dense(x, units=hidden_sizes[-1], activation=output_activation)

    def loss_function(inputs):
        xadv = tf.identity(inputs) # Such that dG/d_xadv = dG/d_x
        G = mlp(xadv, hidden_sizes=(16,1))
        #G, learning_rate = model(xadv)
        dG = tf.gradients(G, xadv)[0]
        # d2G,_ = tf.hessians(G, xadv) # Regularizer
        print('xadv, G and dG shape')
        print(xadv.shape)
        print(G.shape)
        print(dG.shape)
        print((G-xadv[:,:,:,0:1]).shape)
        print('done')

        loss_sum = tf.reduce_mean(((sigma*((G-xadv[:,:,:,0:1])*dG[:,:,:,0:1])) +
            (((xadv[:,:,:,0:1]*G)-beta*xadv[:,:,:,1:2])*dG[:,:,:,1:2]) +
                                   G +
                                   (xadv[:,:,:,0:1]*(xadv[:,:,:,1:2]-gamma))
                                    # + tf.math.scalar_mul(epsilon,(d2Gx+d2Gz))  # Regularizer
                                   )**2)
        print(loss_sum.shape)
        return loss_sum, G # , learning_rate

    """ Optimization setup
    """
    # saver = tf.train.Saver()

    # loss_sum,G,learning_rate = loss_function(x_ph)
    loss_sum,G = loss_function(x_ph)
    # optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate, beta1=0.9,beta2=0.999,epsilon=1e-08)
    optimizer = tf.train.AdamOptimizer(lr, beta1=0.9,beta2=0.999,epsilon=1e-08)
    grads_and_vars = optimizer.compute_gradients(loss_sum)
    op = optimizer.apply_gradients(grads_and_vars,global_step=tf.Variable(0, trainable=False))

    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())


    for epoch in range(epochs): # Minimize loss_sum til |loss-prev_loss|<tol

        db = dataset(reuse = False, **args)
        inds = np.arange(db.shape[0])
        np.random.shuffle(inds)

        for index, ind in enumerate(inds):
            # Sample a batch from cartesian envelope of domain
            # X_sample = sample_domain(batch_size)
            # Feed X to x_ph, do optimization
            #result = sess.run([op, loss_sum, G], feed_dict={x_ph:X_sample, learning_rate:lr(epoch)})
            # result = sess.run([op, loss_sum, G], feed_dict={x_ph:X_sample})
            result = sess.run([op, loss_sum, G], feed_dict={x_ph:db[ind]})
            # if index % 100 == 0:
            #     print("Index = {:5d}; Residuals={: 5.10E};".format(index, result[1]))
            
        if epoch % 1 == 0:
            print("Epoch = {:5d}; Residuals={: 5.10E};".format(epoch, result[1]))
            # Save model weights
            # ...
        if np.abs(result[1]) < 1e-9:
            cprint("Converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
            print("Residual={: 5.10E}.".format(result[1], nrm, err))
            return result[-1]

    cprint("Failed to converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')

    # saver.save(sess, save_path=_SAVE_PATH, global_step=_global_step)

if __name__ == '__main__':
    n_samples_ = 1<<10
    batch_size_ = 1<<5
    sub_nx_, sub_nz_ = 16, 16
    patch_x_, patch_z_ = 1., 1.
    xmin_, xmax_ = -10, 10
    zmin_, zmax_ = -10, 10
    db_save_path_ = './'
    epsilon_ = 0.1
    sigma_ = 10
    beta_ = 8./3
    gamma_ = 25
    solve_lorenz_pde(n_samples=n_samples_, batch_size=batch_size_, 
                     xmin=xmin_, xmax=xmax_, zmin=zmin_, zmax=zmax_, 
                     sub_nx=sub_nx_, sub_nz=sub_nz_, patch_x=patch_x_, patch_z=patch_z_,
                     db_save_path=db_save_path_,
                     epsilon=epsilon_, sigma=sigma_, beta=beta_, gamma=gamma_, epochs=1<<128, lr = 1e-3)









