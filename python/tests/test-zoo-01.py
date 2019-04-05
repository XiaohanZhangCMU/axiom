import sys
sys.path.append('../../lib/')
from zoo import Animal

a = Animal("zoo")
print(a.get_name())
