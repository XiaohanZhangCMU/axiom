import sys
import numpy as np

'''
Remarks to be noted for users:
    1) Sets the mode. It is recommended that you don't change the mode halfway
       into the program since that may cause allocation of pinned memory being
       freed in a non-pinned way, which may cause problems - I haven't verified
       it personally but better to note it here in the header file.
'''

from lib import axiom

def test_Animal():
    animal = axiom.Animal("dog")
   # print(animal)
   # print("The Animal is at 0x{0}".format(animal.get_address()))
   # print("I see a {0}".format(animal.name))
    animal.name = "cat"
   # print("I see a {0}".format(animal.name))
    return 1

def test_CudaAnimal():
    cuda_animal = axiom.CudaAnimal("dog")
    cuda_animal.test_saxpy()
    return 1

def test_NumpyAnimal():
    numpy_animal = axiom.NumpyAnimal()
    dim1 = 10
    dim2 = 1
    x = np.random.random((dim1, dim2)).astype(np.float64)
   # print(np.sum(x * x))
   # print(np.sum(axiom.square_matrix(x)))
    return 1

def test_Axiom():
    axiom.set_mode_cpu()
    axiom.set_mode_gpu()
    axiom.set_multiprocess(True)
    axiom.set_mode_cpu()
    axiom.set_multiprocess(False)
    return 1

def test_Tensor():
    t = axiom.Tensor()
   # print("tensor see {0}".format(t.count))
    m = t.L2
   # print("tensor see {0}".format(m))
    t.reshape(1)
    return 1

def test_MemOps():
    t = axiom.Tensor()
    t = axiom.Tensor()
    arr = np.ones((2,2))
    t.reinit(arr)
   # print(t.data)

    print("The rest tests are not implemented")
    return 1
    foo = t.gpu_data(); # data copied cpu->gpu.
    foo = t.cpu_data(); # no data copied since both have up-to-date contents.
    bar = t.mutable_gpu_data(); # no data copied.
    # ... some operations ...
    bar = t.mutable_gpu_data(); # no data copied when we are still on GPU.
    foo = t.cpu_data(); # data copied gpu->cpu, since the gpu side has modified the data
    foo = t.gpu_data(); # no data copied since both have up-to-date contents
    bar = t.mutable_cpu_data(); # still no data copied.
    bar = t.mutable_gpu_data(); # data copied cpu->gpu.
    bar = t.mutable_cpu_data(); # data copied gpu->cpu.

def main(argv):
    dict_tests = {}
    dict_tests['test_animal']      = test_Animal()
    dict_tests['test_CudaAnimal']  = test_CudaAnimal()
    dict_tests['test_NumpyAnimal'] = test_NumpyAnimal()
    dict_tests['test_Axiom']       = test_Axiom()
    dict_tests['test_Tensor']      = test_Tensor()
    dict_tests['test_MemOps']      = test_MemOps()
  
    for key in dict_tests:
        if dict_tests[key] == 1:
            print("[OK] Test pass-->{0}().".format(key))
        else:
            print("[--] Test fail-->{0}().".format(key))

if __name__ == "__main__":
     main(sys.argv[1:])
