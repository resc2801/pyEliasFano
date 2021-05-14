import math
from typing import List

from bitarray import bitarray
from bitarray.util import ba2int, int2ba
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
        assert all(numbers[i] <= numbers[i + 1] for i in range(len(numbers) - 1)), ValueError("Input list not sorted!")

        self.max_id = (len(numbers) - 1)
        self._size = math.ceil(math.log2(max(numbers) / len(numbers)))
        superiors, inferiors = list(zip(*[divmod(n, 1 << self._size) for n in sorted(numbers)]))

        self._inferiors = bitarray()
        for n in inferiors:
            self._inferiors += bitarray(str(bin(n))[2:].zfill(self._size))

        self._superiors = bitarray()
        for i, n in enumerate(superiors):
            self._superiors += bitarray(inverted_unary(n - superiors[i - 1] if i else n))

    def select(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure.
        :param k: index of integer to be reconstructed.
        :return: k-th stored integer
        """
        assert 0 <= k <= self.max_id, IndexError("Index out of range.")

        inferior = self._get_inferior(k)

        # To compute the higher part we need to perform a select_1(i) - i
        superior = self._superiors.search((int2ba(1)))[k] - k

        return (superior << self._size) + inferior

    def nextGEQ(self, x):
        """
        Return the smallest integer stored in this Elias-Fano structure that is greater or equal than x.
        :param x: integer
        :return: min{y ∈ EF : y ≥ x}
        """
        assert x <= self.select(self.max_id), ValueError("Given value not in this EF structure.")

        if x == self.select(self.max_id):
            return self.select(self.max_id)
        elif x <= self.select(0):
            return self.select(0)
        else:
            sup, inf = divmod(x, 1 << self._size)
            p = int2ba(sup).search((int2ba(0)))[0] - sup

            #        if p < 0:
            #            logging.warn("∀y ∈ EF : x > y")
            #            return self.select(self.max_id)

            for i in range(p, self.max_id + 1):
                if self.select(i) >= x:
                    return self.select(i)

    def nextLEQ(self, x: int):
        """
        Return the largest integer stored in this Elias-Fano structure that is smaller or equal than x.
        :param x: integer
        :return: max{y ∈ EF : x ≥ y}
        """
        assert x >= self.select(0), ValueError("Given value not in this EF structure.")

        if x >= self.select(self.max_id):
            return self.select(self.max_id)
        else:
            return self.select(max(0, self.rank(self.nextGEQ(x)) - 1))

    def rank(self, x: int) -> int:
        """
        Return the index within the Elias-Fano structure for given integer x.
        :param x: integer
        :return: index of x to be reconstructed.
        """
        assert self.select(0) <= x <= self.select(self.max_id), ValueError("Given value not in this EF structure.")
        # TODO: implement rank using bitarray!
        k = -1
        sup, inf = divmod(x, 1 << self._size)
        for c in self._superiors:
            if c == 0:
                sup -= 1
            else:
                k += 1
                if sup == 0 and inf == self._get_inferior(k):
                    return k

        raise ValueError(f"Given element {x} does not exist in structure.")

    def _get_inferior(self, k: int) -> int:
        """
        Return the inferior part of the k-th integer.
        :param k: index
        :return: inferior part of k-th integer stored in this Elias-Fano structure
        """
        assert 0 <= k <= self.max_id, IndexError("Index out of range.")

        inferior = self._inferiors[(self._size * k):(self._size * (k + 1))]
        if len(inferior) == 0:
            return 0
        else:
            return ba2int(inferior)


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

