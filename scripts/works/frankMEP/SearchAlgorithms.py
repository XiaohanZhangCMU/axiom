""" Search Algorihtms: 1) Greedy search. 2) Deep Q-learning.
"""
from MDobj import MDobj
from utility import *
from View import Viewer
from MDobj import MDobj

class GreedySearch(object):
    def __init__(self, cohesv, strain, dirname):
        self.cohesv = cohesv
        self.strain = strain
        self.dirname = dirname
        self.MEP = {}

    #Greedy search for each new atom to add on: always choose the lowest energy state.
    def search(self, swobj, stateA, stateB):
        step0 = saveinter = 0
        nucleus = swobj.reset()
        self.MEP[bits2str(nucleus2bits(nucleus,swobj.totIdx))] = swobj.eval(nucleus)

        print("search I am here 0")

        # Generate a series of ellipse loops, relax with strain
        while True:
            bdyatoms = find_nbr_atoms(swobj.nbrlist, nucleus, swobj.totIdx)
            #assert(bdyatoms.size %2 == 0), "bdyatoms.size = {0}".format(bdyatoms.size)
            #print("bdyatoms = {0}".format(bdyatoms))

            # Find the atom that has lowest energy if added to nucleus
            MINenergy = 1e8
            MINatom_I = -1
            MINatom_J = -1
            bdyatoms_u = np.intersect1d(swobj.pairs[:,0], bdyatoms)
            bdyatoms_d = np.intersect1d(swobj.pairs[:,1], bdyatoms)

            # Make sure loop through the plane having more atoms around nucleus
            nneighbors = bdyatoms_u if len(bdyatoms_u)>len(bdyatoms_d) else bdyatoms_d

            nucleus0 = np.copy(nucleus) # save nucleus before Iter

            if saveinter == 36:
                if 0:
                    tmparr = np.array(bdyatoms).astype(int)
                    red =   [1.0, 0.0, 0.0, 1.0]
                    green = [0.0, 1.0, 0.0, 1.0]
                    blue =  [0.0, 0.0, 1.0, 1.0]
                    plotlist = np.concatenate((nucleus0, tmparr))
                    colorlist = np.vstack((np.tile(red,(len(nucleus0),1)), np.tile(blue,(len(tmparr),1))))
                    view = Viewer(swobj, 300, 300, plotlist, colorlist)
                    view.rendering()
                    self.sw.sleep()

            for atom_I in nneighbors:
                if atom_I not in stateB:
                    continue;
                i,j = np.where(swobj.pairs==atom_I)
                print(swobj.pairs.shape)
                #assert i == 0 or j == 0, "atom_I is either in u or d bdyatoms"
                #if j == 0:
                #    atom_J = swobj.pairs[i[0],1]
                #if i == 0:
                #    atom_J = swobj.pairs[j[0],0]
                atom_J = swobj.pairs[i[0],1-j[0]]

                nucleus_sub_I = np.copy(nucleus)
                nucleus = np.append(nucleus, [atom_I, atom_J])
                nucleus, energy, done, info = swobj.step(nucleus)

                print("search I am here 1")
                nucleus = np.copy(nucleus_sub_I) # save nucleus before subIter

                if energy < MINenergy:
                    MINenergy = energy
                    MINatom_I = atom_I
                    MINatom_J = atom_J

                with open(swobj.dirname+"B.log", "a") as fp:
                    fp.write("nucleus size = {0}; (atom I, atom J) = ({1}); energy = {2}\n".format(len(nucleus), (atom_I, atom_J), energy))
                    #fp.write(str(nucleus))
                    #fp.write('\n')

            if (nucleus2bits(nucleus,swobj.totIdx)==nucleus2bits(stateB,swobj.totIdx)).all():
                print("Reach stateB successfully. GreedySearch() returns.")
                break
            else:
                assert(MINatom_I >=0 and MINatom_J >=0)

            with open(swobj.dirname+"B.log", "a") as fp:
                fp.write("------------------------------\n")
                fp.write("nucleus size = {0}; (atom I, atom J) = ({1}); MINenergy = {2}\n".format(len(nucleus), (MINatom_I, MINatom_J), MINenergy))
                #fp.write(str(nucleus))
                #fp.write('\n')
                fp.write("------------------------------\n")

            nucleus = np.append(nucleus0, [MINatom_I,MINatom_J])

            with open(swobj.dirname+"potential.dat", "a") as fp:
                fp.write(str(MINenergy)+"\n")
            if 1:
                self.save_path_node(swobj, nucleus, saveinter)
                saveinter+=1

            self.MEP[bits2str(nucleus2bits(nucleus,swobj.totIdx))] = MINenergy

            save_obj(self.MEP, swobj.dirname+'mep')

            step0 += 1 # end of for step in range(Nmax)

        return self.MEP

    #def visualize_path():
    #    (X,Y,Z) = order_path(self.MEP)
    #    plt.plot(X,Y,'*-')
    #    plt.draw()

    def write_path(path, name):
        (N,Y,X) = order_path(self.MEP)
        with open(name, "w") as fp:
            for val in N:
                fp.write(str(val)+" ")
            fp.write("\n")
            for val in Y:
                fp.write(str(val)+" ")
            fp.write("\n")

    def order_path():
        X = np.array(list(self.MEP.keys()))
        Y = np.array(list(self.MEP.values()))
        N = [ ]
        for item in X:
            N.append(nbits1(str2bits(item)))
        N = np.array(N)
        inds = N.argsort()
        Y = Y[inds]
        N = N[inds]
        X = X[inds]
        return (N, Y, X)

    def visualize_path_node(self, swobj, nucleus):
        red =   [1.0, 0.0, 0.0, 1.0]
        green = [0.0, 1.0, 0.0, 1.0]
        blue =  [0.0, 0.0, 1.0, 1.0]
        plotlist =np.extract(np.abs(swobj.SR[:,2])>0.375, np.arange(swobj.sw.NP))
        atomlist = np.concatenate((plotlist,nucleus))
        colorlist = np.vstack((np.tile(red,(len(plotlist),1)), np.tile(blue,(len(nucleus),1))))
        view = Viewer(swobj, 300, 300, atomlist, colorlist)
        view.rendering()

    def save_path_node(self, swobj, nucleus, saveinter):
        swobj.sw.finalcnfile=swobj.dirname + "/img_"+str(saveinter)+".cfg"
        swobj.sw.writeatomeyecfg(swobj.sw.finalcnfile)

        swobj.sw.finalcnfile=swobj.dirname + "/pathimg_"+str(saveinter)+".cfg"
        swobj.sw.freeallatoms()
        swobj.fixed.fill(1)
        swobj.fixed[nucleus] = 0
        print("I am here 1")
        swobj.sw.removefixedatoms()
        swobj.sw.writeatomeyecfg(swobj.sw.finalcnfile)

        swobj.fixed.fill(0)
        print("I am here 2")
        swobj.sw.NP = swobj.NP0
        swobj.sw.SR1toSR()
        print("I am here 3")
        swobj.sw.refreshnnlist()
        print("I am here 4")


class DQNSearch(object):
    def __init__(self, swobj, cohesv, strain, dirname,
                n_actions,
                n_features,
                learning_rate=0.01,
                reward_decay=0.9,
                e_greedy=0.9,
                replace_target_iter=300,
                memory_size=500,
                batch_size=32,
                e_greedy_increment=None,
                output_graph=False,):
        self.cohesv = cohesv
        self.strain = strain
        self.dirname = dirname
        self.n_actions = n_actions
        self.n_features = n_features
        self.lr = learning_rate
        self.gamma = reward_decay
        self.epsilon_max =e_greedy
        self.replace_target_iter = replace_target_iter
        self.memory_size = memory_size
        self.batch_size = batch_size
        self.epsilon_increment = e_greedy_increment
        self.epsilon = 0 if e_greedy_increment is not None else self.epsilon_max

        # total learning step
        self.learn_ste_counter = 0

        # Initialize zero memory [s, a, r, s_]
        self.memory = np.zeros((self.memory_size, n_features*2+2))
        # consist of [target_net, evaluate_net]
        self._build_net()
        t_params = tf.get_collection('target_net_params')
        e_params = tf.get_collection('eval_net_params')
        self.replace_target_op = [tf.assign(t,e) for t, e in zip(t_params, e_params)]

        self.sess = tf.Session()

        if output_graph:
            tf.summary.FileWriter("logs/", self.sess.graph)

        self.sess.run(tf.global_variables_initializer())
        self.cost_hist = []

    def _build_net(self):
        self.s = tf.placeholder(tf.float32, [None, self.n_features], name='s') # input
        self.q_target = tf.placeholder(tf.float32, [None, self.n_actions], name='Q_target') # loss calculating

        with tf.variable_scope('eval_net'):
            # collection names
            c_name, n_l1, w_initializer, b_initializer = \
                    ['eval_net_params', tf.GraphKeys.GLOBAL_VARIABLES], 10, \
                    tf.random_normal_initializer(0., 0.3), tf.constant_initializer(0.1) # config of layers

            # first layer
            with tf.variable_scope('l1'):
                w1 = tf.get_variable('w1', [self.n_features, n_l1], initializer=w_initializer, collections=c_names)
                b1 = tf.get_variable('b1', [1, n_l1], initializer=b_initializer, collections=c_names)
                l1 = tf.nn.relu(tf.matmul(self.s, w1)+b1)
            # second layer
            with tf.variable_scope('l2'):
                w2 = tf.get_variable('w2', [n_l1, self.n_actions], initializer=w_initializer, collections=c_names)
                b2 = tf.get_variable('b2', [1, self.n_actions], initializer=b_initializder, collections=c_names)
                self.q_next = tf.matmul(l1, w2)+b2

            with tf.variable_scope('loss'):
                self.loss = tf.reduce_mean(tf.squared_difference(self.q_target, self.q_eval))
            with tf.variable_scope('train'):
                self._train_op = tf.train.RMSPropOptimizer(self.lr).minimizer(self.loss)

        self.s_ = tf.placeholder(tf.float32, [None, self.n_features], name='s_') # input
        with tf.varaible_scope('target_net'):
            c_names = ['target_net_params', tf.GraphKeys.GLOBAL_VARIABLES]

            # first layer
            with tf.varialbe_scope('l1'):
                w1 = tf.get_variable('w1', [self.n_features, n_l1], initializer=w_initializer, collections=c_names)
                b1 = tf.get_variable('b1', [1, n_l1], initializer=b_initializer, collections=c_names)
                l1 = tf.nn.relue(tf.matmul(self.s_,w1)+b1)

            # second layer
            with tf.varialbe_scope('l2'):
                w2 = tf.get_variable('w2', [n_l1, self.n_actions], initializer=w_initializer, collections=c_names)
                b2 = tf.get_variable('b2', [1, self.n_actions], initializer=b_initializer, collections=c_names)
                self.q_next = tf.matmul(l1, w2) + b2


    def store_transition(self, s, a, r, s_):
        if not hasattr(self, 'memory_counter'):
            self.memory_counter = 0

        transition = np.hstack((s, [a, r], s_))

        # replace the old memory with new memory
        index = self.memory_counter % self.memory_size
        self.memory[index, :] = transition

        self.memory_counter += 1

    def choose_action(self, observation):
        # to have batch dimension when feed into tf placeholder
        observation = observation[np.newaxis,:]

        if np.random.uniform() < self.epsilon:
            # forward feed observation and get q value
            actions_value = self.sess.run(self.q_eval, feed_dict={self.s : observation})
            action = np.argmax(actions_value)
        else:
            action = np.random.randint(0, self.n_actions)
        return action

    def learn(self):
        # check to replace target parameters
        if self.learn_step_counter % self.replace_target_iter == 0:
            self.sess.run(self.replace_target_op)
            print('\ntarget_params_replaced')

        # sample batch memory from all memory
        if self.memory_counter > self.memory_size:
            sample_index = np.random.choice(self.memory_size, size = self.batch_size)
        else:
            sample_index = np.random.choice(self.memory_counter, size = self.batch_size)
        batch_memory = self.memory[sample_index, :]

        q_next, q_eval = self.sess.run(
                [self.q_next, self.q_eval],
                feed_dict = {
                    self.s_:batch_memory[:, -self.n_features:], # fixed params
                    self.s :batch_memory[:, :self.n_features],  # newest params
                    })
        # change q_target w.r.t q_eval`s action
        q_target = q_eval.copy()

        batch_index = np.arange(self.batch_size, dtype=np.int32)
        eval_act_index = batch_memory[:, self.n_features].astype(int)
        reward = batch_memory[:, self.n_features+1]

        q_target[batch_index, eval_act_index] = reward + self.gamma * np.max(q_next, axis=1)

        # train eval network
        _, self.cost = self.sess.run([self._train_op, self.loss], feed_dict={self.s: batch_memory[:, :self.n_features], self.q_target: q_target})

        self.cost_his.append(self.cost)

        # increasing epsilon
        self.epsilon = self.epsilon + self.epsilon_increment if self.epsilon < self.epsilon_max else self.epsilon_max
        self.learn_step_counter += 1

    def search(self, Nmax):
        step = 0

        for episode in range(Nmax):
            observation = env.reset()

            while True:
                env.render()

                # RL choose action based on observation
                action = self.choose_action(observation)

                # RL take action and get next observation and reward
                observation_, reward, done = env.step(action)

                self.store_transition(observation, action, reward, observation_)

                if (step > 200) and (step %5 == 0):
                    self.learn()

                observation = observation_

                if done:
                    break
                step +=1

        print('Search finished')

    def plot_cost(self):
        import matplotlib.pyplot as plt
        plt.plot(np.arange(len(self.cost_his)), self.cost_his)
        plt.ylabel('Cost')
        plt.xlabel('training steps')

