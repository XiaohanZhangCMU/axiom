import numpy as np
import sys
sys.path.append('../../lib/')
import tensor

'''
Remarks to be noted for users:
    1) Sets the mode. It is recommended that you don't change the mode halfway
       into the program since that may cause allocation of pinned memory being
       freed in a non-pinned way, which may cause problems - I haven't verified
       it personally but better to note it here in the header file.
'''

def test_Tensor():
  error = 0
  t = tensor.Tensor()
# print("tensor see {0}".format(t.count))
  m = t.L2
# print("tensor see {0}".format(m))
  t.Reshape(1)
  return error

def test_MemOps():
  error = 0
  t1 = tensor.Tensor()
  t2 = tensor.Tensor()
  arr = np.ones((2,2))
  t1.reinit(arr)
  t2.reinit(arr)

# k = tensor.Tensor()
# k.assign( t1[0] )

  if sum(t1.data - t2.data).sum() != 0:
    error += 1
  arr = np.ones((50,1))
  t1.reinit(arr)
  t3 = t1
  if sum(t1.data - t3.data).sum() != 0:
    error += 1
  arr = np.ones(t1.shape)
  if sum(t1.data - arr).sum() != 0:
    error += 1
  t1 = t2
  if sum(t3.data - arr).sum() != 0:
    error += 1
  return error

def test_GPU():
  error = 0
  t = tensor.Tensor()
  #foo = t.gpu_data(); # data copied cpu->gpu.
  #foo = t.cpu_data(); # no data copied since both have up-to-date contents.
  #bar = t.mutable_gpu_data(); # no data copied.
  # ... some operations ...
  #bar = t.mutable_gpu_data(); # no data copied when we are still on GPU.
  #foo = t.cpu_data(); # data copied gpu->cpu, since the gpu side has modified the data
  #foo = t.gpu_data(); # no data copied since both have up-to-date contents
  #bar = t.mutable_cpu_data(); # still no data copied.
  #bar = t.mutable_gpu_data(); # data copied cpu->gpu.
  #bar = t.mutable_cpu_data(); # data copied gpu->cpu.
  return 0


def main(argv):
  dict_tests = {}
  dict_tests['test_Tensor']      = test_Tensor()
  dict_tests['test_MemOps']      = test_MemOps()
  dict_tests['test_GPU']         = test_GPU()

  for key in dict_tests:
    if dict_tests[key] == 0:
      print("[OK] Test pass-->{0}().".format(key))
    else:
      print("[--] Test fail-->{0}().".format(key))

if __name__ == "__main__":
  main(sys.argv[1:])

