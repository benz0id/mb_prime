import timeit
from scan_for_binding_regions import main

time = timeit.timeit(main, number=1)

print('Runtime:', time, 'seconds')
