import math
from functools import reduce
from typing import List

import numpy as np
from bitarray import bitarray
from bitarray.util import ba2int
from bitarray.util import int2ba, zeros, count_n


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
        :param numbers: list of integers SORTED IN ASCENDING ORDER
        """

        # sequence length
        self._n = len(numbers)

        # size of universe
        self._u = 2**(math.ceil(math.log2(numbers[-1])))

        # number of upper bits per sequence element
        self._upper_bits = math.ceil(math.log2(self._n))

        # number of lower bits per sequence element
        self._lower_bits = math.floor(math.log2(self._u / self._n))

        # upper bits of each sequence element in negated unary encoding
        self._superiors = self._encode_upper_bits(numbers)

        # lower bits of each sequence element in fixed width representation
        self._inferiors = self._encode_lower_bits(numbers)

    def select(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure.
        Note that we use 1-based indexing!
        :param k: index of integer to be reconstructed.
        :return: k-th stored integer
        """
        if not(0 < k <= self._n):
            raise IndexError(f"Use any k ∈ [1,..,{self._n}].")

        # for lower part simply jump to the corresponding bits in self._inferiors
        inferior = ba2int(self._inferiors[(self._lower_bits * (k-1)):(self._lower_bits * k)])

        # To compute the higher part we need to perform a select_1(k) - k on self._superiors
        superior = count_n(self._superiors, k) - k    # returns lowest index i for which a[:i].count() == n

        return (superior << self._lower_bits) | inferior

    def rank(self, x: int) -> int:
        """
        Return the index within the Elias-Fano structure for given integer x.
        :param x: integer
        :return: index of x to be reconstructed.
        """
        if not(self.select(1) <= x <= self.select(self._n)):
            raise ValueError(f"{x} ∉ EF.")

        # split x into upper_bits and lower_bits
        sup, inf = ((x >> self._lower_bits) & ((2 ** self._upper_bits) - 1), (x & ((2 ** self._lower_bits) - 1)))

        # fetch the correct bucket of lower bits
        b = self._superiors[:list(self._superiors.itersearch(zeros(1)))[sup]].count(1)
        if sup > 0:
            a = self._superiors[:list(self._superiors.itersearch(zeros(1)))[sup-1]].count(1)
        else:
            a = 0

        # loop through the bucket and search for lower bits of x
        for k in range(a, b):
            if ba2int(self._inferiors[(self._lower_bits * k):(self._lower_bits * (k+1))]) == inf:
                return k+1

        raise ValueError(f"{x} ∉ EF.")

    def nextGEQ(self, x):
        """
        Return the smallest integer stored in this Elias-Fano structure that is greater or equal than x.
        :param x: integer
        :return: min{y ∈ EF : y ≥ x}
        """
        if not(x <= self.select(self._n)):
            raise ValueError(f"∄y: min{{y ∈ EF: y ≥ {x}}}.")

        if x <= self.select(1):
            return self.select(1)

        elif x == self.select(self._n):
            return self.select(self._n)

        elif self.select(1) < x < self.select(self._n):
            # split x into upper_bits and lower_bits
            sup, inf = ((x >> self._lower_bits) & ((2 ** self._upper_bits) - 1), (x & ((2 ** self._lower_bits) - 1)))

            # fetch the correct bucket of lower bits
            b = self._superiors[:list(self._superiors.itersearch(zeros(1)))[sup]].count(1)
            if sup > 0:
                a = self._superiors[:list(self._superiors.itersearch(zeros(1)))[sup - 1]].count(1)
            else:
                a = 0

            # loop through the bucket and search for lower_bits >= x
            for k in range(a, b):
                if ba2int(self._inferiors[(self._lower_bits * k):(self._lower_bits * (k + 1))]) >= inf:
                    return self.select(k + 1)

            # if bucket did not contain nextGEQ(x), we take select(b+1)
            return self.select(b+1)

    def nextLEQ(self, x: int):
        """
        Return the largest integer stored in this Elias-Fano structure that is smaller or equal than x.
        :param x: integer
        :return: max{y ∈ EF : x ≥ y}
        """
        if not(x >= self.select(1)):
            raise ValueError(f"∄y: max{{y ∈ EF : {x} ≥ y}}")

        if x >= self.select(self._n):
            return self.select(self._n)
        else:
            return x if x == self.nextGEQ(x) else self.select(max(1, self.rank(self.nextGEQ(x)) - 1))

    def _encode_upper_bits(self, numbers: List[int]):
        """
        Encodes the upper bits of each sequence element into a bitarray
        :param numbers: list of integers SORTED IN ASCENDING ORDER
        :return: bitarray containing all upper bits in negated unary encoding
        """
        # the upper bits for each sequence element
        superiors = [((v >> self._lower_bits) & ((2 ** self._upper_bits) - 1)) for v in numbers]

        # represent upper bits for each sequence element in negated unary
        (bins, counts) = np.unique(superiors, return_counts=True)
        negated_unary_upper_bits = [0] * (2 ** self._upper_bits)
        for i in range(0, bins.shape[0]):
            negated_unary_upper_bits[bins[i]] = (((2 ** counts[i]) - 1) << 1).item()

        # return negated unary upper bits for all sequence items as single bitarray
        return reduce(lambda a, b: a + b, list(map(int2ba, negated_unary_upper_bits)), bitarray())

    def _encode_lower_bits(self, numbers: List[int]):

        # the lower bits for each sequence element
        inferiors = [(v & ((2 ** self._lower_bits) - 1)) for v in numbers]

        # store lower bits in fixed width into bitarray
        fixed_width_lower_bits = [ zeros(self._lower_bits - v.bit_length()) + int2ba(v) for v in inferiors]

        return reduce(lambda a, b: a + b, fixed_width_lower_bits, bitarray())


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
