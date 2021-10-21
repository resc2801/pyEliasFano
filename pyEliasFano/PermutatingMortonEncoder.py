from itertools import chain, zip_longest, accumulate
from typing import List
from functools import reduce
from sympy.combinatorics import Permutation
from more_itertools import take


class PermutatingMortonEncoder:
    """
    This Encoder accepts n-dimensional integer vectors, applies a given permutation to each input dimension and
    subsequently shuffles the permutated input components into a Morton code.
    """

    def __init__(self, permutation_per_dim: List[Permutation] = None):
        """
        Initializes the PermutatingMortonEncoder with given component permutations.
        Identity permutations allow standard Morton encoding.
        """
        self._permutation_per_dim = permutation_per_dim

        # this permutation will shuffle the concatentation of the components' binary representations into a Morton code
        # Virtuous Death!
        self._shuffle = Permutation(
            list(
                filter(lambda x: x is not None,
                       chain.from_iterable(
                           zip_longest(*map(lambda start_stop: range(start_stop[1] - start_stop[0], start_stop[1]),
                                            zip(list(map(lambda x: x.size, permutation_per_dim)),
                                                accumulate(list(map(lambda x: x.size, permutation_per_dim))))
                                            )
                                       )
                       )
                       )
            )
        )

    def encode(self, data_point: List[int]) -> int:
        """
        Encodes a n-dimensional integer vector into a Morton code.
        """
        # apply permutation per data_point component
        permuted_data_point = reduce(
            lambda a, b: a + b,
            ["".join(
                self._permutation_per_dim[i]((("{0:0%db}" % self._permutation_per_dim[i].size).format(data_point[i]))))
                for i in range(len(self._permutation_per_dim))]
        )

        # shuffle the components' bits into Morton code
        return int(reduce(lambda a, b: a + b, self._shuffle(permuted_data_point)), 2)

    def decode(self, morton_code: int) -> List[int]:
        """
        Decodes a given Morton code into its respective n-dimensional integer vector.
        """
        # get the bits per components
        bits_per_component = list(map(lambda x: x.size, self._permutation_per_dim))

        # de-shuffle the Morton bits
        bit_iter = iter((~self._shuffle)(("{0:0%db}" % sum(bits_per_component)).format(morton_code)))

        # apply inverse permutation per component
        return [int("".join((~(self._permutation_per_dim[i]))(perm)), 2) for i, perm in
                enumerate(
                    map(lambda no_bits: "".join(list(take(no_bits, bit_iter))),
                        bits_per_component
                        )
                )]
