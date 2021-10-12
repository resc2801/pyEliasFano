import math
from collections import Counter
from itertools import accumulate, islice, chain
from typing import List, Iterator

from more_itertools import locate, first, nth, unzip


class EliasFano:
    """
    An Elias-Fano structure represents a monotone non-decreasing sequence of n integers from the universe [0 . . . m)
    occupying 2n+n⌈log2(m/n)⌉ bits.

    It supports:
     - select(k): nearly constant time access to the k-th element,
     - rank(x): fast access to the index within the structure for given integer x.
     - nextGEQ(x): fast access to the smallest integer of the sequence that is greater or equal than a given x
     - nextLEQ(x): fast access to the largest integer of the sequence that is smaller or equal than a given x
    """

    def __init__(self, sorted_integers: List[int]):
        """
        Construct an Elias-Fano structure for the given sorted list of integers.
        :param sorted_integers: list of integers SORTED IN ASCENDING ORDER
        """
        # sequence length
        self._n = len(sorted_integers)

        # size of universe
        self._u = 2 ** max(1, sorted_integers[-1].bit_length())

        # number of upper bits per sequence element
        self._upper_bits = math.ceil(math.log2(self._n))

        # number of lower bits per sequence element
        self._lower_bits = max(0, math.floor(math.log2(self._u / self._n)))

        # upper bits of each sequence element in negated unary encoding
        # lower bits of each sequence element in fixed width representation
        self._encode(sorted_integers)

    def select(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure.
        :param k: index of integer to be reconstructed.
        :return: k-th stored integer
        """
        if not (0 <= k < self._n):
            raise IndexError(f"Use any k ∈ [0,..,{self._n - 1}].")

        # for lower part simply take self._inferiors[k], defaults to 0 in case self._lower_bits == 0
        inferior = nth(iter(self._inferiors), k, 0)

        # the higher part is index of the bucket with accumulated popCount >= k, defaults to 0 iff self._upper_bits == 0
        superior = first(locate(islice(iter(self._superiors_prefixSums), 1, None), lambda cnt: cnt > k), 0)

        return (superior << self._lower_bits) | inferior

    def rank(self, x: int) -> List[int]:
        """
        Return the indices within the Elias-Fano structure for given integer x.
        :param x: integer
        :return: list of indices of x in the Elias-Fano structure.
        """

        # split x into upper_bits and lower_bits
        sup_x = (x >> self._lower_bits) & ((2 ** self._upper_bits) - 1)
        inf_x = (x & ((2 ** self._lower_bits) - 1))

        # count elements in EF structure with sup_x as upper_bits
        if self._upper_bits > 0:
            a, b = (self._superiors_prefixSums[sup_x], self._superiors_prefixSums[sup_x + 1])
        else:
            a, b = (0, len(self._inferiors))

        if self._lower_bits > 0:
            return list(map(lambda k: a + k,
                            locate(islice(self._inferiors, a, b),
                                   pred=lambda inf: inf == inf_x)))
        else:
            return list(range(a, b))

    def nextGEQ(self, x) -> int:
        """
        Return the smallest integer stored in this Elias-Fano structure that is greater or equal than x.
        :param x: integer
        :return: min{y ∈ EF : y ≥ x}
        """
        if not (x <= self.select(self._n - 1)):
            raise ValueError(f"∄y: min{{y ∈ EF: y ≥ {x}}}.")

        if x <= self.select(0):
            return self.select(0)

        # split x into upper_bits and lower_bits
        sup_x = (x >> self._lower_bits) & ((2 ** self._upper_bits) - 1)
        inf_x = (x & ((2 ** self._lower_bits) - 1))

        # count elements in EF structure with sup_x as upper_bits
        a = max(0, self._superiors_prefixSums[sup_x])
        b = min(self._superiors_prefixSums[sup_x + 1], len(self._inferiors))

        k = min(first(map(lambda i: a + i, # loop through bucket and search for inf_bits >= x
                          locate(islice(self._inferiors, a, b),
                                 pred=lambda inf: inf >= inf_x)), self._n),
                b)  # if bucket did not contain nextGEQ(x), we take select(b)

        return self.select(k)

    def nextLEQ(self, x: int) -> int:
        """
        Return the largest integer stored in this Elias-Fano structure that is smaller or equal than x.
        :param x: integer
        :return: max{y ∈ EF : x ≥ y}
        """
        if not (x >= self.select(0)):
            raise ValueError(f"∄y: max{{y ∈ EF : {x} ≥ y}}")

        if x >= self.select(self._n - 1):
            return self.select(self._n - 1)

        if x == self.nextGEQ(x):
            return x

        return self.select(max(0, first(self.rank(self.nextGEQ(x))) - 1))

    def bit_length(self):
        """
        Number of bits needed.
        """
        return sum(map(lambda ones: ones + 1, self._superiors)) + self._lower_bits * len(self._lower_bits)

    def compression_ratio(self):
        """
        Compression ratio is defined as the ratio between the uncompressed size and compressed size
        """
        return (self._n * (math.log2(self._u))) / self.bit_length()

    def _encode(self, sorted_integers: List[int]):
        """
        Compresses a monotone non-decreasing integers lists by using Elias-Fano encoding.
        :param sorted_integers: list of integers SORTED IN ASCENDING ORDER
        """
        inferiors_iter, superiors_iter = unzip(
            map(lambda x: ((x & ((2 ** self._lower_bits) - 1)),
                           ((x >> self._lower_bits) & ((2 ** self._upper_bits) - 1))),
                iter(sorted_integers)))

        if self._lower_bits > 0:
            # the lower bits for each sequence element
            self._inferiors = list(inferiors_iter)
        else:
            # empty list if we do not use lower_bits
            self._inferiors = []

        if self._upper_bits > 0:
            # init all upper_bit buckets with 0
            self._superiors = [0] * (2 ** self._upper_bits)

            # encode the count of `superior` bits (in negated unary representation when serializing)
            for (superior, count) in Counter(superiors_iter).items():
                self._superiors[superior] = count
        else:
            # empty list if we do not use upper_bits
            self._superiors = []

        # auxiliary prefix_sum array to speed up rank computation
        self._superiors_prefixSums = [0] + list(accumulate(self._superiors))

    def __getitem__(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure using subscript operator.
        """
        return self.select(k)

    def __len__(self) -> int:
        """
        Return number of integers stored in the Elias-Fano structure.
        """
        return self._n

    def __iter__(self) -> Iterator[int]:
        """
        Support for __iter__ and next
        """
        # iterate elements in _superiors, _superiors[k] is the numbers of elements in _inferiors
        # to fetch and combine with k as their upper_half
        _inferiors_iter = iter(self._inferiors)
        return chain.from_iterable(
            map(lambda sup: [(sup[0] << self._lower_bits) + inf for inf in islice(_inferiors_iter, sup[1])],
                enumerate(self._superiors)))


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
