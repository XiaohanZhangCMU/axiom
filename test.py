from axiomlib import add
add(2, 3)

from axiomlib import Pet
my_dog = Pet('Pluto', 5)
my_dog.get_name()

print(my_dog.get_hunger())

my_dog.go_for_a_walk()
print(my_dog.get_hunger())
