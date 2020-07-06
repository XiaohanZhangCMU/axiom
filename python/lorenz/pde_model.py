import tensorflow as tf
import numpy as np

class Lorenz_Model:
    def __init__(self, **config):
        sigma = config['sigma']
        beta = config['beta']
        gamma = config['gamma']
        Bbar = config['Bbar']

        self.x_ph = x_ph = tf.placeholder(dtype=np.float64,
                shape=(None, 2), name='input')
        xadv = tf.identity(x_ph) # Such that dG/d_xadv = dG/d_x
        self.G = G = self.mlp_policy(xadv, hidden_sizes=(8,2,1))
        dG = tf.gradients(G, xadv)[0]
        # d2G,_ = tf.hessians(G, xadv) # Regularizer
        # self.loss_sum = loss_sum = tf.reduce_sum(((sigma*((G-xadv[:,0:1])*dG[:,0:1])) + (((xadv[:,0:1]*G)-beta*xadv[:,1:2])*dG[:,1:2]) + G + (xadv[:,0:1]*(xadv[:,1:2]-gamma)))**2)
        # Scaled form
        self.loss_sum = loss_sum = tf.reduce_sum(((sigma*((G-xadv[:,0:1])*dG[:,0:1])) + ((Bbar*xadv[:,0:1]*G-beta*xadv[:,1:2])*dG[:,1:2]) + G + (xadv[:,0:1]*(Bbar*xadv[:,1:2]-gamma)))**2)
        # + tf.math.scalar_mul(epsilon,(d2Gx+d2Gz))  # Regularizer

        print('xadv, G and dG shape')
        print(xadv.shape)
        print(G.shape)
        print(dG.shape)
        print((G-xadv[:,0:1]).shape)
        print(loss_sum.shape)
        print('done')

        self.learning_rate = tf.placeholder(tf.float32, shape=[], name='learning_rate')
        optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate,beta1=0.9,beta2=0.999,epsilon=1e-08)
        grads_and_vars = optimizer.compute_gradients(loss_sum)
        self.op = optimizer.apply_gradients(grads_and_vars,global_step=tf.Variable(0, trainable=False))

    def mlp_policy(self, x, hidden_sizes=(32,1), activation=tf.nn.sigmoid, output_activation=tf.nn.sigmoid):
        for h in hidden_sizes[:-1]:
            x = tf.layers.dense(x, units=h, activation=activation, kernel_initializer=tf.contrib.layers.xavier_initializer())
        return tf.layers.dense(x, units=hidden_sizes[-1], activation=output_activation, kernel_initializer=tf.contrib.layers.xavier_initializer())

    def cnn_policy(self, x_image, sub_nx, feat_dim=1):
        _IMAGE_SIZE = sub_nx
        _IMAGE_CHANNELS = 2
        with tf.name_scope('main_params'):
            learning_rate = tf.placeholder(tf.float32, shape=[], name='learning_rate')
        with tf.variable_scope('conv1') as scope:
            conv = tf.layers.conv2d(inputs=x_image,filters=32,kernel_size=[3, 3],
                padding='SAME',activation=tf.nn.relu)
            conv = tf.layers.conv2d(inputs=conv,filters=64,kernel_size=[3, 3],
                padding='SAME',activation=tf.nn.relu)
            #pool = tf.layers.max_pooling2d(conv, pool_size=[2, 2], strides=2, padding='SAME')
            #drop = tf.layers.dropout(pool, rate=0.25, name=scope.name)
            drop = conv
        with tf.variable_scope('conv2') as scope:
            conv = tf.layers.conv2d(inputs=drop,filters=128, kernel_size=[3, 3],
                padding='SAME',activation=tf.nn.relu)
            #pool = tf.layers.max_pooling2d(conv, pool_size=[2, 2], strides=2, padding='SAME')
            pool = conv
            conv = tf.layers.conv2d(inputs=pool,filters=128,kernel_size=[2, 2],
                padding='SAME',activation=tf.nn.relu)
            #pool = tf.layers.max_pooling2d(conv, pool_size=[2, 2], strides=2, padding='SAME')
            #drop = tf.layers.dropout(pool, rate=0.25, name=scope.name)
            drop = conv
        with tf.variable_scope('fully_connected') as scope:
            #flat = tf.reshape(drop, [-1, 4 * 4 * 128])
            fc = tf.layers.dense(inputs=drop, units=feat_dim, activation=tf.nn.relu)
        return fc, learning_rate

