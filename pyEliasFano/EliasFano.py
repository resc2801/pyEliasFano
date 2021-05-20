import math
from typing import List
from gmpy2 import xmpz, mpz, t_divmod_2exp, t_div_2exp, t_mod_2exp
from itertools import compress, islice
import operator
import gc
from bitarray import bitarray
from bitarray.util import ba2int
from unary_coding import inverted_unary


class EliasFano:
    """
    An Elias-Fano structure represents a monotone non-decreasing sequence of n integers from the universe [0 . . . m)
    occupying 2n+n⌈log2(m/n)⌉ bits.

    It supports:
     - select(k): nearly constant time access to the k-th element,
     - rank(x): access to the index within the structure for given integer x.
     - nextGEQ(x): fast access to the smallest integer of the sequence that is greater or equal than a given x
     - nextLEQ(x): fast access to the largest integer of the sequence that is smaller or equal than a given x
    """

    def __init__(self, numbers: List[int]):
        """
        Construct an Elias-Fano structure for the given sorted list of integers.
        :param numbers: sorted list of integers
        """

        self._min_val = numbers[0]
        self._max_val = numbers[-1]
        self._m = max(numbers)
        self._n = len(numbers)
        self._upper_bits = math.ceil(math.log2(self._n))
        self._lower_bits = math.ceil(math.log2(self._m / self._n))
        self._inferiors = xmpz()  # n*log2(m/n) bits
        self._superiors = xmpz()  # 2n bit

        sups, infs = list(zip(*[t_divmod_2exp(n, self._lower_bits) for n in numbers]))

        for i, n in enumerate(infs):
            self._inferiors[i * self._lower_bits:(i + 1) * self._lower_bits] = n[0:self._lower_bits]

        superiors = bitarray()
        for i, n in enumerate(sups):
            superiors += bitarray(inverted_unary(int(n) - sups[i - 1] if i else int(n)))
        superiors.reverse()
        self._superiors = xmpz(ba2int(superiors))

    def select(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure.
        :param k: index of integer to be reconstructed.
        :return: k-th stored integer
        """
        assert 0 <= k < self._n, IndexError("Index out of range.")

        # for lower part simply jump to the corresponding bits in self._inferiors
        inferior = int(self._inferiors[(self._lower_bits * k):(self._lower_bits * (k + 1))])

        # To compute the higher part we need to perform a select_1(k) - k on self._superiors
        superior = next(islice(compress(range(0, 2 * self._n), list(self._superiors.iter_bits())), k, k + 1)) - k

        return (superior << self._lower_bits) + inferior

    def rank(self, x: int) -> int:
        """
        Return the index within the Elias-Fano structure for given integer x.
        :param x: integer
        :return: index of x to be reconstructed.
        """
        assert self._min_val <= x <= self._max_val, ValueError("Value not in this EF structure.")

        k = -1
        sup, inf = t_divmod_2exp(x, self._lower_bits)

        for c in self._superiors.iter_bits():
            if c:
                k += 1
                if sup == 0 and inf == int(self._inferiors[(self._lower_bits * k):(self._lower_bits * (k + 1))]):
                    return k
            else:
                sup -= 1

        raise ValueError(f"{x} ∉ EF.")

    def nextGEQ(self, x):
        """
        Return the smallest integer stored in this Elias-Fano structure that is greater or equal than x.
        :param x: integer
        :return: min{y ∈ EF : y ≥ x}
        """
        assert x <= self._max_val, ValueError(f"∄y: min{{y ∈ EF: y ≥ {x}}}.")

        if x <= self.select(0):
            return self._min_val

        elif x == self.select(self._n - 1):
            return self._max_val

        elif self.select(0) < x < self.select(self._n - 1):
            sup = t_div_2exp(x, self._lower_bits)
            p = next(islice(compress(range(0, 2 * self._n), map(operator.not_, self._superiors.iter_bits())), sup-1, sup)) - (sup-1)
            for i in range(p, self._n):
                if self.select(i) >= x:
                    return self.select(i)

    def nextLEQ(self, x: int):
        """
        Return the largest integer stored in this Elias-Fano structure that is smaller or equal than x.
        :param x: integer
        :return: max{y ∈ EF : x ≥ y}
        """
        assert x >= self._min_val, ValueError(f"∄y: max{{y ∈ EF : {x} ≥ y}}")

        if x >= self._max_val:
            return self._max_val
        else:
            return x if x == self.nextGEQ(x) else self.select(max(0, self.rank(self.nextGEQ(x)) - 1))


def load(file_path: str):
    """
    Load Elias-Fano structure from disk.
    :param file_path: file storing an Elias-Fano structure
    :return: an Elias-Fano structure
    """
    import os
    from pickle import load

    assert os.path.exists(file_path) and os.path.isfile(file_path), IOError("File path invalid or does not exist.")

    obj = load(open(file_path, "rb"))
    assert isinstance(obj, EliasFano), ValueError("Given file does not store an Elias-Fano structure.")

    return obj


def save(ef_structure: EliasFano, file_path: str):
    """
    Save given Elias-Fano structure to disk.
    :param ef_structure: an Elias-Fano structure
    :param file_path: file to hold the EF structure
    :return: None
    """
    from pickle import dump

    assert isinstance(ef_structure, EliasFano), ValueError("Given object is not an Elias-Fano structure.")

    dump(ef_structure, open(file_path, "wb"))
