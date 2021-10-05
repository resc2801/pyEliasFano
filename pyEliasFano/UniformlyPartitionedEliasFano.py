import logging
import math
from itertools import chain
from typing import List, Iterator
from pyEliasFano import EliasFano

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


class UniformlyPartitionedEliasFano:
    """
    The basic idea is to partition the sequence S into m/b chunks of b consecutive integers each,
    except possibly the last one.

    The first level is an Elias-Fano representation of the sequence L obtained by juxtaposing the last element of each
    chunk (i.e., L = S[0], S[b], S[2b], ... ).

    The second level is the collection of the chunks of S, each represented with Elias-Fano.
    """

    def __init__(self, numbers: List[int], b: int):

        # sequence length
        self._n = len(numbers)

        # size of universe
        self._u = 2 ** (numbers[-1].bit_length())

        # chunk size, except possibly the last one
        self._b = b

        #  partition the sequence S into m/b chunks of b consecutive integers each, except possibly the last one
        chunks = [numbers[i:i + self._b] for i in range(0, len(numbers), self._b)]

        # 1st level is an EF representation of the sequence L obtained by juxtaposing the 1st element of each chunk
        L = [chunk[0] for chunk in chunks]

        # elements of the jth chunk can be rewritten in a smaller universe by subtracting L[j] from each element
        chunks = [list(map(lambda x: x - L[j], chunk)) for (j, chunk) in enumerate(chunks)]

        logging.info('Constructing 1st level of EF index for %s chunks.' % len(chunks))
        self._level_1 = EliasFano(L)

        logging.info('Constructing 2nd level of EF index for %s chunks.' % len(chunks))
        self._level_2 = {key: EliasFano(chunk) for (key, chunk) in enumerate(chunks, start=0)}

        logging.info('Done.')

    def select(self, i: int) -> int:
        """
        Return i-th integer stored in this uniformly-partitioned EF structure. Indexing is zero-based.
        :param i: index of integer to be reconstructed.
        :return: i-th integer in sequence
        """
        if not(0 <= i < self.__len__()):
            raise IndexError(f"Use any i ∈ [0,..,{self.__len__() - 1}].")

        # Let j be the index of the chunk containing the ith element of S (i.e., j = floor(i/b))
        # Let k be its offset within this chunk (i.e., k = i mod b)
        (j, k) = (math.floor(i / self._b), i % self._b)

        # Knowing L[j] and b suffices for accessing the kth element of the jth chunk.
        e = self._level_2[j][k]

        # If e is the value at this position, then we conclude that the value S[i] is equal to L[j]+e
        return self._level_1[j] + e

    def rank(self, x: int) -> int:
        """
        Return the index i within the uniformly-partitioned EF structure for given integer x.
        :param x: integer
        :return: index of x to be reconstructed.
        """
        if not(self.select(0) <= x <= self.select(self._n - 1)):
            raise ValueError(f"{x} ∉ EF.")

        try:
            # let j be the index of the chunk that might contain x
            j = self._level_1.rank(self._level_1.nextLEQ(x))

            # the j-th chunk contains e = x - L[j] iff x is contained within the uniformly-partitioned EF structure
            e = x-self._level_1.select(j)

            # let k be the index of e = x-L[j] within the j-th chunk
            k = self._level_2[j].rank(e)

            return (j * self._b) + k

        except ValueError:
            raise ValueError(f"{x} ∉ EF.")

    def bit_length(self):
        """
        Return the combined bit lengths of the first and second level EF structures
        of this uniformly-partitioned EF structure.
        """
        return self._level_1.bit_length() + sum(ef2.bit_length() for ef2 in (self._level_2.values()))

    def compression_ratio(self):
        """
        Return the data compression ratio achieved with this uniformly-partitioned EF structure.
        Data compression ratio is defined as the ratio between the uncompressed size and compressed size.
        """
        return (self._n * (math.log2(self._u))) / self.bit_length()

    def __getitem__(self, i: int) -> int:
        self.select(i)

    def __iter__(self) -> Iterator[int]:
        return chain.from_iterable(map(lambda t: [t[0]+v for v in t[1]],
                                       zip(iter(self._level_1),
                                           iter(self._level_2.values()))))

    def __len__(self) -> int:
        return sum(len(ef2) for ef2 in (self._level_2.values()))


