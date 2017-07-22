from lib import axiom

#assert 'hello' in dir(axiom)
#assert callable(axiom.hello)
#print(axiom.hello())

animal = axiom.Animal("dog")

print(animal)
print("The C++ object is at 0x{0}".format(animal.get_address()))
print("I see a {0}".format(animal.name))

animal.name = "cat"
print("I see a {0}".format(animal.name))

cuda_animal = axiom.CudaAnimal("dog")
cuda_animal.test_saxpy()
