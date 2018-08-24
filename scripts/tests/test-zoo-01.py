import sys
sys.path.append('../../lib/')
from axiomlib import Animal

a = Animal("zoo")
print(a.get_name())
