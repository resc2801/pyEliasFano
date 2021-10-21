import math
from collections import Counter
from itertools import accumulate, islice, chain, starmap
from typing import List, Iterator, Tuple

from more_itertools import locate, first, nth, unzip

# (value_bits, dont_care_bits)
BitPattern = Tuple[int, int]


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

        # size of universe (NOTE: python assumes int(0).bit_length() == 0)
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
        Returns the sequence element 'x' such that 'k' sequence predecessors are smaller or equal to 'x'.
        :param k: index of integer to be reconstructed.
        :return: k-th stored integer
        """
        if not (0 <= k < self._n):
            raise IndexError(f"Index %s ∉ [0,..,{self._n - 1}]." % k)

        # for lower part simply take self._inferiors[k], defaults to 0 in case self._lower_bits == 0
        inferior = nth(iter(self._inferiors), k, 0)

        # the higher part is index of the bucket with accumulated popCount >= k, defaults to 0 iff self._upper_bits == 0
        superior = first(locate(iter(self._superiors_prefixSums), lambda cnt: cnt > k), 0)

        return (superior << self._lower_bits) | inferior

    def rank(self, x: int) -> List[int]:
        """
        Return the indices within the Elias-Fano structure for given integer x.
        :param x: integer
        :return: list of indices of x in the Elias-Fano structure.
        """

        # split x into upper_bits and lower_bits
        sup_x, inf_x = self._split(x)

        if self._upper_bits > 0:
            inferior_range = range(self._superiors_prefixSums[sup_x] - self._superiors[sup_x],
                                   self._superiors_prefixSums[sup_x])
            if self._lower_bits > 0:
                return list(filter(lambda i: self._inferiors[i] == inf_x, inferior_range))
            else:
                return list(inferior_range)
        else:
            if self._lower_bits > 0:
                return locate(iter(self._inferiors), lambda inf: inf == inf_x)

        raise ValueError("%s ∉ EF structure." % x)

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
        sup_x, inf_x = self._split(x)

        inferior_range = range(self._superiors_prefixSums[sup_x] - self._superiors[sup_x],
                               self._superiors_prefixSums[sup_x])

        k = min(first(filter(lambda i: self._inferiors[i] >= inf_x, inferior_range), self._n),
                self._superiors_prefixSums[sup_x])

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

    def match(self, value, ignore):

        # split value_bits into upper_bits and lower_bits
        sup_value, inf_value = self._split(value)

        # split dont_care_bits into upper_bits and lower_bits
        sup_ignore, inf_ignore = self._split(ignore)

        # filter matching upper halves in self._superiors
        if self._upper_bits > 0:
            sup_matches = list(filter(lambda idx: (idx & sup_ignore) == (sup_value & sup_ignore),
                                 locate(self._superiors, pred=lambda cnt: cnt > 0)))

            # for each matching upper half, we pair it with its matching lower halves
            matching_elements = list(map(lambda sup_match: (sup_match,
                                                       filter(
                                                           lambda inf: (inf & inf_ignore) == (inf_value & inf_ignore),
                                                           islice(self._inferiors,
                                                                  self._superiors_prefixSums[sup_match]-self._superiors[sup_match],
                                                                  self._superiors_prefixSums[sup_match]))),
                                    sup_matches))

            # reconstruct stored integer
            return chain.from_iterable(
                map(lambda pair: map(lambda inf: (pair[0] << self._lower_bits) + inf, pair[1]), matching_elements)
            )
        else:
            return filter(lambda low: (low & inf_ignore) == (inf_value & inf_ignore), self._inferiors)

    def bit_length(self):
        """
        Number of bits needed.
        """
        return sum(map(lambda ones: ones + 1, self._superiors)) + self._lower_bits * len(self._inferiors)

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

        # encode the lower bits for each sequence element; [] if we do not use lower_bits
        self._inferiors = []
        if self._lower_bits > 0:
            self._inferiors = list(inferiors_iter)

        # encode the upper bits for each sequence element; [] if we do not use upper_bits
        self._superiors = []
        if self._upper_bits > 0:
            # init all upper_bit buckets with 0
            self._superiors = [0] * (2 ** self._upper_bits)

            # encode the count of `superior` bits (in negated unary representation when serializing)
            for (superior, count) in Counter(superiors_iter).items():
                self._superiors[superior] = count

        self._superiors_prefixSums = list(accumulate(self._superiors))

    def _split(self, x: int) -> Tuple[int, int]:
        """
        Split binary representation of integer x into upper and lower half.
        :param x: integer
        :return: tuple with upper and lower half of x
        """
        return (x >> self._lower_bits) & ((2 ** self._upper_bits) - 1), (x & ((2 ** self._lower_bits) - 1))

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
        if self._upper_bits:
            if self._lower_bits:
                # iterate elements in _superiors, _superiors[k] is the numbers of elements in _inferiors
                # to fetch and combine with k as their upper_half
                _inferiors_iter = iter(self._inferiors)
                return chain.from_iterable(
                    map(lambda idx: [(idx << self._lower_bits) + inf for inf in
                                     islice(_inferiors_iter, self._superiors[idx])],
                        locate(self._superiors, pred=lambda cnt: cnt > 0)))
            else:
                return map(lambda sup: (sup << self._lower_bits), locate(self._superiors, pred=lambda cnt: cnt > 0))
        else:
            if self._lower_bits:
                return iter(self._inferiors)
            else:
                raise ValueError("Empty index!")


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
