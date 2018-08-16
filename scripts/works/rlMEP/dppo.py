""" A practice script of OpenAI`s Proximal Policy Optimization (PPO) with multithreading
    https://arxiv.org/abs/1707.06347
    
    Distributing workers in parallel to collect data, then stop worker`s roll-out and train PPO on collected data. 
    Restart workers once PPO is updated. 

    The global PPO updating rule is adopted from DeepMind`s paper (DPPO):
    Emergence of Locomotion Behaviors in Rich Environments (Google Deepmind): https://arxiv.org/abs/1707.02286
    Origin: https://morvanzhou.github.io/tutorials
"""

import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import gym, threading, queue

EP_MAX = 1000        # max experience
EP_LEN = 200
N_WORKERS = 4        # num of threads or workers in parallel 
GAMMA = 0.9          # reward discount factor
A_LR = 0.0001        # learning rate for actor
C_LR = 0.0002        # learning rate for critic
MIN_BATCH_SIZE = 64  # minimum batch size for PPO updating 
UPDATE_STEP = 10     # loop update operation n-steps
EPSILON = 0.2        # clipping surrogate objective
GAME = 'Pendulum-v0' # OpenAI gym regisry identifiers
S_DIM, A_DIM = 3, 1  # state and action dimension 


class PPO(object):
    def __init__(self):
        self.sess = tf.Session()
        self.tfs = tf.placeholder(tf.float32, [None, S_DIM], 'state') # None means the axis does not have restriction, i.e. input feature size x batch size

        #critic
        l1 = tf.layers.dense(self.tfs, 100, tf.nn.relu) # input size, output size
        self.v = tf.layers.dense(l1, 1)
        self.tfdc_r = tf.placeholder(tf.load32, [ None, 1], 'discounted_r')
        self.advantage = self.tfdc_r - self.v
        self.closs = tf.reduce_mean(tf.square(self.advantage))
        self.ctrain_op = tf.train.AdamOptimizer(C_LR).minimize(self.closs)
        
        # actor
        pi, pi_params = self._build_anet('pi', trainable=True)
        oldpi, oldpi_params = self._build_anet('oldpi', trainable=False)
        self.sample_op = tf.squeeze(pi.sample(1), axis=0) #choosing actions
        self.update_oldpi_op = [ oldp.assign(p) for p, oldp in zip(pi_params, oldpi_params)]
        self.tfa = tf.placeholder(tf.float32, [None, A_DIM], 'action')
        self.tfadv = tf.placeholder(tf.float32, [None, 1], 'advantage')
        ratio = pi.prob(self.tfa) / (oldpi.prob(self.tfa) + 1e-5)
        surr = ratio * self.tfadv

        self.aloss = -tf.reduce_mean(tf.minimum(  # clipped surrogate objective
            surr, 
            tf.clip_by_value(ratio, 1. - EPSILON, 1. + EPSILON)) * self.tfadv)

        self.atrain_op = tf.train.AdamOptimizer(A_LR).minimize(self.aloss)
        self.sess.run(tf.global_variables_initializer())



































