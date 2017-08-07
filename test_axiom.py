from lib import axiom
import numpy as np

'''
Remarks to be noted for users:

    Sets the mode. It is recommended that you don't change the mode halfway
    into the program since that may cause allocation of pinned memory being
    freed in a non-pinned way, which may cause problems - I haven't verified
    it personally but better to note it here in the header file.

Todo: 
    compile .....
'''

'''
Test class Animal
'''
animal = axiom.Animal("dog")
print(animal)
print("The Animal is at 0x{0}".format(animal.get_address()))
print("I see a {0}".format(animal.name))
animal.name = "cat"
print("I see a {0}".format(animal.name))

'''
Test class CudaAnimal
'''
cuda_animal = axiom.CudaAnimal("dog")
cuda_animal.test_saxpy()

'''
Test class NumpyAnimal
'''
numpy_animal = axiom.NumpyAnimal()
dim1 = 10
dim2 = 1
x = np.random.random((dim1, dim2)).astype(np.float64)
print(np.sum(x * x))
print(np.sum(axiom.square_matrix(x)))

'''
Test class Axiom
'''
#axiom.set_mode_cpu()
