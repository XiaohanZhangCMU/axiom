import os
import numpy as np
import tensorflow as tf
import functools
import collections
import multiprocessing

def load_params(params_txt):
    def isFloat(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    f = open(params_txt,'r')
    d = {}
    for line in f:
        if '=' not in line:
            continue
        s = line.split('=')[1]
        if '\\n' in line:
            s = s[:-2]
        s = s.strip()
        if isFloat(s):
            if s.isnumeric(): # integer
                s = int(s)
            else:
                s = float(s)
        else:
            s = s[1:-1] # remove ' ' or " " around string
        d[line.split('=')[0].strip() ] = s
    f.close()
    return d

def get_session(config=None):
    """Get default session or create one with a given config"""
    sess = tf.get_default_session()
    if sess is None:
        sess = make_session(config=config, make_default=True)
    return sess

def make_session(config=None, num_cpu=None, make_default=False, graph=None):
    """Returns a session that will use <num_cpu> CPU's only"""
    if num_cpu is None:
        num_cpu = int(os.getenv('RCALL_NUM_CPU', multiprocessing.cpu_count()))
    if config is None:
        config = tf.ConfigProto(
            allow_soft_placement=True,
            inter_op_parallelism_threads=num_cpu,
            intra_op_parallelism_threads=num_cpu)
        config.gpu_options.allow_growth = True

    if make_default:
        return tf.InteractiveSession(config=config, graph=graph)
    else:
        return tf.Session(config=config, graph=graph)

def save_variables(save_path, variables=None, sess=None):
    import joblib
    sess = sess or get_session()
    variables = variables or tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES)

    ps = sess.run(variables)
    save_dict = {v.name: value for v, value in zip(variables, ps)}
    dirname = os.path.dirname(save_path)
    if any(dirname):
        os.makedirs(dirname, exist_ok=True)
    joblib.dump(save_dict, save_path)

def load_variables(load_path, variables=None, sess=None):
    import joblib
    sess = sess or get_session()
    variables = variables or tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES)

    loaded_params = joblib.load(os.path.expanduser(load_path))
    restores = []
    if isinstance(loaded_params, list):
        assert len(loaded_params) == len(variables), 'number of variables loaded mismatches len(variables)'
        for d, v in zip(loaded_params, variables):
            restores.append(v.assign(d))
    else:
        for v in variables:
            restores.append(v.assign(loaded_params[v.name]))

    sess.run(restores)
