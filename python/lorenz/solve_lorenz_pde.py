""" Use tensorflow to solve PDE of lorenz-G equation
    1) Formulate least square loss
    2) Sample PDE domain patchwise
"""
import sys
sys.path.append('../../lib/')
from termcolor import colored, cprint
import lorenz
import tensorflow as tf

def fc(x, scope, nh, *, init_scale=1.0, init_bias=0.0):
    with tf.variable_scope(scope):
        nin = x.get_shape()[1].value
        w = tf.get_variable("w", [nin, nh], initializer=ortho_init(init_scale))
        b = tf.get_variable("b", [nh], initializer=tf.constant_initializer(init_bias))
        return tf.matmul(x, w)+b

def mlp(num_layers=2, num_hidden=64, activation=tf.tanh, layer_norm=False):
    def network_fn(X):
        h = tf.layers.flatten(X)
        for i in range(num_layers):
            h = fc(h, 'mlp_fc{}'.format(i), nh=num_hidden, init_scale=np.sqrt(2))
            h=activation(h)
        return h
    return network_fn

def solve_lorenz_pde(xmin, xmax, zmin, zmax, epsilon):
    sub_nx, sub_nz = 10, 10
    
    x_ph = tf.placeholder(shape=(None, sub_nx, sub_nz,2),name='input')
    model = mlp(num_layers=2,num_hidden=32)

    def sample_domain(batch=32):
        xs = uniform(xmin, xmax, batch)
        zs = uniform(zmin, zmax, batch)
        xn = linspace(xs, xs+patch_size_x, sub_nx)
        zn = linspace(zs, zs+patch_size_z, sub_nz)
        return np.concat(xn,zn,axis=1)

    def loss_function(inputs):
        xadv = tf.identity(inputs) # Such that dG/d_xadv = dG/d_x
        G = model.netwrok_fn(xadv)
        dG,_ = tf.gradients(G, xs=xadv)
        # d2G,_ = tf.hessians(G, xs=xadv) # Regularizer
        loss_sum = tf.reduce_mean(tf.math.square( 
            	tf.math.scalar_mul(sigma, tf.math.multiply((G-xadv), Gx)) 
            	+ tf.math.multiply((tf.math.scalar_mul(x,G)-tf.math.scalar_mul(beta,z)),Gz)  
                + G  
                + tf.math.multiply(x, (z-gamma)) 
                # + tf.math.scalar_mul(epsilon,(d2Gx+d2Gz))  # Regularizer
                ))
        return loss_sum

    """ Optimization setup
    """
    lr = 1e-3  # learning rate
    tol = 1e-9 # convergence tol
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
            start = batch * batch_size
            end = min(n_sample, start+batch_size) 
            # Feed X to x_ph, do optimization
            result = sess.run([op, loss_sum, G], feed_dict={x_ph:X_data[start:end]})

        if epoch % 100 == 0:
            print("Epoch = {:5d}; Residuals={: 5.10E};".format(epoch, result[1]))
            # Save model weights 
            # ...
        if err < 1e-9:
            cprint("Converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
            print("Residual={: 5.10E}.".format(result[1], nrm, err))
            return result[-1]

    cprint("Failed to converged in {:d} epochs!!!".format(epoch), 'green', 'on_red')
