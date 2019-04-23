"""
    Use tensorflow to solve PDE of lorenz-G equation
    1) Formulate least square loss
    2) Sample PDE domain patchwise

    Questions to resolve:
    1) Run time shape dimension does not match?
    2) Normalize input patch to have values [-1,1] for each channel?
    3) Add laplacian regularizer to PDE.
"""
import sys
from termcolor import colored, cprint
import numpy as np
import tensorflow as tf


def solve_lorenz_pde(xmin, xmax, zmin, zmax, epsilon, sigma, beta, gamma):
    sub_nx, sub_nz = 16, 16
    patch_x, patch_z = 1., 1.

    x_ph = tf.placeholder(dtype=np.float64, shape=(None, sub_nx, sub_nz, 2), name='input')
    model = mlp(num_layers=2, num_hidden=32)

    def mlp(x, hidden_sizes=(32,), activation=tf.tanh, output_activation=None):
        for h in hidden_sizes[:-1]:
            x = tf.layers.dense(x, units=h, activation=activation)
        return tf.layers.dense(x, units=hidden_sizes[-1], activation=output_activation)

    def sample_domain(batch_size=32):
        xs = np.random.uniform(xmin, xmax, batch_size)
        zs = np.random.uniform(zmin, zmax, batch_size)
        return np.array([np.stack(np.meshgrid(np.linspace(xi-0.5*patch_x, xi+0.5*patch_x, sub_nx),
            np.linspace(zi-0.5*patch_z, zi+0.5*patch_z, sub_nz)), axis=-1)
            for (xi,zi) in zip(xs,zs)])

    def loss_function(inputs):
        xadv = tf.identity(inputs) # Such that dG/d_xadv = dG/d_x
        G = mlp(xadv, hidden_sizes=[32,1])
        dG = tf.gradients(G, xadv)[0]
        # d2G,_ = tf.hessians(G, xadv) # Regularizer
        print('G and dG shape')
        print(xadv.shape)
        print(G.shape)
        print(dG.shape)

        loss_sum = tf.reduce_mean(((sigma*((G-xadv[:,0])*dG[:,0])) +
                                   (((xadv[:,0]*G)-beta*xadv[:,1])*dG[:,1]) +
                                   G +
                                   (xadv[:,0]*(xadv[:,1]-gamma))
                                    # + tf.math.scalar_mul(epsilon,(d2Gx+d2Gz))  # Regularizer
                                   )**2)
        print(loss_sum.shape)
        return loss_sum, G

    """ Optimization setup
    """
    lr = 1e-3  # learning rate
    tol = 1e-9 # convergence tol
    batch_size = 1<<5
    epochs = 1<<20

    optimizer = tf.train.AdamOptimizer(lr)
    grads_and_vars = optimizer.compute_gradients(loss_sum)
    op = optimizer.apply_gradients(grads_and_vars,global_step=tf.Variable(0, trainable=False))

    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    sess.run(tf.local_variables_initializer())

    for epoch in range(epochs): # Minimize loss_sum til |loss-prev_loss|<tol

        # Sample a batch from cartesian envelope of domain
        ind = np.arange(n_sample)
        np.random.shuffle(ind)

        for batch in range(n_batch):
            X_sample = sample_domain(batch_size)
            # Feed X to x_ph, do optimization
            batch_loss = loss_function(x_ph)
            result = sess.run([op, batch_loss, G], feed_dict={x_ph:X_sample})

        if epoch % 100 == 0:
            print("Epoch = {:5d}; Residuals={: 5.10E};".format(epoch, result[1]))
            # Save model weights
            # ...
        if err < 1e-9:
            cprint("Converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
            print("Residual={: 5.10E}.".format(result[1], nrm, err))
            return result[-1]

    cprint("Failed to converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')

if __name__ == '__main__':
    solve_lorenz_pde(xmin=-10, xmax=10, zmin=-10, zmax=10, epsilon=0.1)



