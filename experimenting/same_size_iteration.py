from time import time
from typing import List


def smart_incr(lst: List[int], max_val: int) -> bool:
    def getmin() -> int:
        for i in range(len(lst) - 1):
            if lst[i] > 0 and lst[i + 1] < max_val:
                return i
        return -1

    min_ind = getmin()

    if min_ind == -1:
        return False

    lst[min_ind] -= 1
    lst[min_ind + 1] += 1

    sums = sum(lst[:min_ind + 1])
    for i in range(min_ind + 1):
        lst[i] = 0

    i = 0
    while sums != 0:
        amt = min(sums, max_val)
        lst[i] = amt
        sums -= amt
        i += 1

    return True




def foo(lst, max_val) -> float:
    t1 = time()
    num = 1
    inds = [i for i in range(len(lst))]

    def getmin() -> int:
        for i in range(len(lst) - 1):
            if lst[i] > 0 and lst[i + 1] < max_val:
                return i
        print(time() - t1)
        print(lst)
        print(num)
        return -1

    print(lst, '\n')

    while smart_incr(lst, max_val):
        num += 1
    print(time() - t1)
    print(lst)
    print(num)


if __name__ == '__main__':
    # foo([12, 12, 12, 12, 0, 0, 0, 0], 12)
    foo([3, 0, 0], 12)