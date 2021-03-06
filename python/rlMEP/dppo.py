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

    def update(self):
        global GLOBAL_UPDATE_COUNTER
        while not COORD.should_stop():
            if GLOBAL_EP < EP_MAX:
                UPDATE_EVENT.wait()   # wait until get batch of data
                self.sess.run(self.update_oldpi_op) #copy pi to old_pi
                data = [ q.get() for _ in range(q.qsize()) ] #collect data from workers
                data = np.vstack(data)

                s,a,r = data[:,:S_DIM], data[:,S_DIM:S_DIM+A_DIM], data[:,-1:]
                adv = self.sess.run(self.advantage, { self.tfs:s, self.tfa:a, self.tfadv:adv})

                #update actor and critic in a update loop
                [self.sess.run(self.atrain_op, { self.tfs: s , self.tfa: a, self.tfadv: adv}) for _ in range(UPDATE_STEP)]

                [self.sess.run(self.ctrain_op, { slef.tfs: s, self.tfdc_r: r}) for _in range(UPDATE_STEP)]

    def _build_anet(self, name, trainable):
        with tf.variable_scope(name):
            l1 = tf.layers.dense(self.tfs, 200, tf.nn.relu, trainable=trainable)
            mu = 2 * tf.layers.dense(l1, A_DIM, tf.nn.tanh, trainable=trainable)
            sigma = tf.layers.dense(l1, A_DIM, tf.nn.softplus, trainable=trainable)
            norm_dist = tf.distributions.Normal(loc=mu, scale= sigma)
        params = tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES, scope=name)
        return norm_dist, params

    def choose_action(self, s):
        s = s[np.newaxis, :]  # add another dimension to s
        a = self.sess.run(self.sample_op, {self.tfs: s})[0]
        return np.clip(a, -2, 2) # values <-2 will become -2, and those >2 will become 2

    def get_v(self, s):
        if s.ndim<2 : s = s[np.newaxis,:]
        return self.sess.run(self.v, {self.tfs:s})[0, 0]

class Worker(object):
    def __init__(self, wid):
        self.wid = wid     # worker id
        self.env = gym.make(GAME).unwrapped
        self.ppo = GLOBAL_PPO

    def work(self):
        global GLOBAL_EP, GLOBAL_RUNNING_R, GLOBAL_UPDATE_COUNTER
        while not COORD.should_stop():
            s = self.env.reset()    # reset starts a new environment
            ep_r = 0
            buffer_s, buffer_a, buffer_r = [], [], []  # buffer for observation, action, reward

            for t in range(EP_LEN):
                if not ROLLING_EVENT.is_set():        # while global PPO is updating
                    ROLLING_EVENT.wait()              # wait until PPO is updated
                    buffer_s, buffer_a, buffer_r = [], [], [] # clear histroy buffer, use new policy to collect data
                a = self.ppo.choose_action(s)
                """ env.step returns four values
                    observation: state of environemnt
                    reward of previous action.
                    done: True if the episode is finished
                    info: diagnostic information for debug
                """
                s_, r, done, _ = self.env.step(a)
                buffer_s.append(s)
                buffer_a.append(a)
                buffer_r.append((r+8)/8)              # normalize reward, find to be useful
                s = s_
                ep_r += r

                GLOBAL_UPDATE_COUNTER += 1            # count to minimum batch size, no need to wait other workers

                if t == EP_LEN - 1 or GLOBAL_UPDATE_COUNTER >= MIN_BATCH_SIZE:
                    v_s = self.ppo.get_v(s_)
                    discounted_r = []
                    for r in buffer_r[::-1]:
                        v_s = r + GAMMA * v_s_
                        discounted_r.append(v_s_)
                    discounted_r.reverse()

                    bs, ba, br = np.vstack(buffer_s), np.vstack(buffer_a), np.array(discounted_r)[:, np.newaxis]
                    buffer_s, buffer_a, buffer_r = [], [], []
                    q.put(np.hstack((bs, ba, br))
                    if GLOBAL_UPDATE_COUNTER >= MIN_BATCH_SIZE:
                        ROLLING_EVENT.clear()    # stop collecting data
                        UPDATE_EVENT.set()       # globalPPO update

                    if GLOBAL_EP >= EP_MAX:      # stop training
                        COORD.request_stop()
                        break

            # record reward changes, plot later
            if len(GLOBAL_RUNNING_R) == 0: GLOBAL_RUNNING_R.append(ep_r)
            else: GLOBAL_RUNNING_R.append(GLOBAL_RUNNING_R[-1]*0.9+ep_r*0.1)
            GLOBAL_EP += 1
            print('{0:.1f}%'.format(GLOBAL_EP/EP_MAX*100), '|W%i' % self.wid, '|Ep_r: %.2f' % ep_r,)


if __name__ == '__main__':
    GLOBAL_PPO = PPO()
    UPDATE_EVENT, ROLLING_EVENT = threading.Event(), threading.Event()
    UPDATE_EVENT.clear()        # not update now
    ROLLING_EVENT.set()         # start to roll out
    workers = [ Worker(wid=i) for i in range(N_WORKER) ] # initialize N_WORKER workers

    GLOBAL_UPDATE_COUNTER, GLOBAL_EP = 0, 0
    GLOBAL_RUNNING_R = []
    COORD = tf.train.Coordinator() # Coordinator() coordinates termination of a set of threads
    q = queue.Queue()           # workers put data in this queue
    threads = []
    for worker in workers:      # worker threads
        t = threading.Thread(target=worker.work, args=())
        t.start()               # training
        threads.append(t)

    # add a PPO updating thread
    threads.append(threading.Thread(target=GLOBAL_PPO.update,))
    threads[-1].start()
    COORD.join(threads)         # wait for all threads in threads to finish


    #plot reward change and test
    plt.plot(np.arange(len(GLOBAL_RUNNING_R)), GLOBAL_RUNNING_R)
    plt.xlable('Episode'); plt.ylable('Moving reward'); plt.ion(); plt.show();

    env = gym.make('Pendulum-v0')
    while True:
        s = env.reset()
        for t in range(300):
            env.render()
            s = env.step(GLOBAL_PPO.choose_action(s))[0]



































