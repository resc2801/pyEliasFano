import math
from functools import reduce
from typing import List, Iterator
from itertools import islice
import numpy as np
from bitarray import bitarray
from bitarray.util import ba2int
from bitarray.util import int2ba, zeros, count_n
from pyEliasFano import EliasFano
from itertools import chain

import logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


class UniformPartitionedEliasFano:
    """
    The basic idea is to partition the sequence S into m/b chunks of b consecutive integers each,
    except possibly the last one.

    The first level is an Elias-Fano representation of the sequence L obtained by juxtaposing the last element of each
    chunk (i.e., L = S[b − 1], S[2b − 1], . . . , S[m − 1]).

    The second level is the collection of the chunks of S, each represented with Elias-Fano.
    """

    def __init__(self, numbers: List[int], b: int):

        # sequence length
        self._n = len(numbers)

        # size of universe
        self._u = 2 ** (numbers[-1].bit_length())

        logging.info('Constructing uniform partitioned EF index for %s integers.' % len(numbers))
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

        logging.info('Done.' % len(chunks))

    def select(self, i: int) -> int:
        """
        Return i-th integer stored in this uniformly partitioned Elias-Fano structure.
        Note that we use 1-based indexing!
        :param i: index of integer to be reconstructed.
        :return: i-th integer in sequence
        """
        if not(0 <= i < self.__len__()):
            raise IndexError(f"Use any k ∈ [0,..,{self.__len__() - 1}].")

        # Let j be the index of the chunk containing the ith element of S (i.e., j = floor(i/b))
        # Let k be its offset within this chunk (i.e., k = i mod b)
        (j, k) = (math.floor(i / self._b), i % self._b)

        # access L[j−1] and L[j] on the first level to compute the size of the universe U_j of the chunk as L[j]−L[j−1]
        # Knowing U_j and b suffices for accessing the kth element of the jth chunk.
        e = self._level_2[j][k]

        # If e is the value at this position, then we conclude that the value S[i] is equal to L[j]+1+e
        return self._level_1[j] + e

    def __getitem__(self, i: int) -> int:
        self.select(i)

    def __iter__(self) -> Iterator[int]:
        return chain.from_iterable(map(lambda t: [t[0]+v for v in t[1]],
                                       zip(iter(self._level_1),
                                           iter(self._level_2.values()))))

    def __len__(self) -> int:
        return len(self._level_1) + sum(len(ef2) for ef2 in (self._level_2.values()))

    def bit_length(self):
        return self._level_1.bit_length() + sum(ef2.bit_length() for ef2 in (self._level_2.values()))

    def compression_ratio(self):
        return self.bit_length() / (self._n * (math.log2(self._u)))

if __name__ == "__main__":

    import numpy as np
    import csv
    import math
    from pyEliasFano import EliasFano, UniformPartitionedEliasFano
    import logging

    logging.info('Loading RDF from disk.')
    output = []
    with open('/Users/resc01-admin/VisualStudio/morton/data/eswc.tsv') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=' ')

        for row in csvreader:
            output.append([int(row[0]), int(row[1]), int(row[2])])

    triples = np.array(output)
    largest_ids = max(triples[:, 0]), max(triples[:, 1]), max(triples[:, 2])
    component_bitsizes = [int(id).bit_length() for id in largest_ids]

    logging.info('Representing %s triples in binary.' % len(triples))
    inputs = sorted(
        [((s << (sum(component_bitsizes[1:3]))) | (p << sum(component_bitsizes[2:3])) | o) for s, p, o in output])

    foo = UniformPartitionedEliasFano(inputs, 10000)
    print(foo.select(0))
    print(foo.select(27000))
    for i in range(len(foo)):
        assert foo.select(i) == inputs[i], f"Failure for index {i}"
