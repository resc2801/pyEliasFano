import logging
import math
from itertools import accumulate, islice, chain
from itertools import groupby
from operator import itemgetter
from typing import List, Iterator
from more_itertools import locate, first, nth, unzip, take
from pyEliasFano import EliasFano


class MultiLevelEliasFano:
    """
    A multi-level Elias-Fano structure represents a monotone non-decreasing sequence of n integers from the universe [0 . . . m).

    It supports:
     - select(k): nearly constant time access to the k-th element,
    """

    def __init__(self, sorted_integers: List[int], depth: int):
        """
        Construct a multi-level Elias Fano structure for the given sorted list of integers and with given depth.
        :param sorted_integers: list of integers SORTED IN ASCENDING ORDER
        :param depth: depth of the multi-level Elias Fano structure
        """

        # sequence length
        self._n = len(sorted_integers)

        # size of universe (NOTE: python assumes int(0).bit_length() == 0)
        self._u = 2 ** max(1, sorted_integers[-1].bit_length())

        # depth
        self._d = depth

        # number of input bits handled on this level
        self._b = divmod(max(1, sorted_integers[-1].bit_length()), self._d)[0]

        # secondary level
        self._level_2 = {}

        # if we are not at the last level
        if (self._d > 1) and (self._n > 1):
            # for each sequence element in fixed width representation, strip-off a prefix of `b` bits
            # as upper halves and group with its respective suffices (aka lower halves)
            chunks = groupby(map(lambda x: ((x >> max(1, sorted_integers[-1].bit_length()) - self._b),
                                            (x & int("1" * (max(1, sorted_integers[-1].bit_length()) - self._b), 2))),
                                 sorted_integers),
                             key=lambda x: x[0])

            L = []  # stores all upper halves for lvl-1 encoding

            for (prefix_val, suffix_grouper) in chunks:
                L.append(prefix_val)
                # for each prefix, encode its respective suffixes in a lvl-2 (MultiLevel)EliasFano
                suffixes = list(map(itemgetter(1), suffix_grouper))
                if len(suffixes) > (2**self._b):
                    self._level_2[prefix_val] = MultiLevelEliasFano(suffixes, self._d - 1)
                else:
                    self._level_2[prefix_val] = EliasFano(suffixes)

            # next store all prefix bits in a lvl-1 Elias Fano structure
            self._level_1 = EliasFano(L)

        # at the last level we simply encode the sequence into a Elias Fano structure
        else:
            self._level_1 = EliasFano(sorted_integers)

    def select(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure.
        :param k: index of integer to be reconstructed.
        :return: k-th stored integer
        """

        # if we are not at the last level
        if bool(self._level_2):

            # determine the 'h'-th lvl-2 bucket_id containing the k-th element's lower half
            # note that 'h' is a list index - NOT a prefix label
            h = first(locate(accumulate(map(lambda e: len(e),
                                            self._level_2.values())),
                             lambda popcnt: popcnt > k))

            # determine number 'l' of elements contained in lvl-2 buckets [0..h)
            l = sum(map(lambda e: len(e), islice(self._level_2.values(), 0, h)))

            # the '(k-l)'-th element from lvl-2 index with prefix_label['h'] must be the 'k'-th element's lower half
            inf = self._level_2[nth(self._level_2.keys(), h)].select(k - l)

            # get the 'k'-th element's upper half from lvl-1 and left_shift appropriately
            sup = self._level_1.select(h) << (int(math.log2(self._u)) - self._b)

            # return the combined upper and lower halves
            return sup + inf

        # if at last level
        else:
            # we don't have a lvl-2 index and therefore return the 'k-th element
            return self._level_1.select(k)

    def match(self, value, ignore):

        suffix_bit_width = (int(math.log2(self._u)) - self._b)
        (sup_value, sup_ignore) = ((value >> suffix_bit_width), (ignore >> suffix_bit_width))

        lvl1_matches = list(filter(lambda x: (x & sup_ignore) == (sup_value & sup_ignore), iter(self._level_1)))

        if bool(self._level_2):
            (inf_value, inf_ignore) = (value & ((2**suffix_bit_width) - 1), ignore & ((2**suffix_bit_width) - 1))
            return chain.from_iterable(
                map(lambda lvl1_match: [(lvl1_match << suffix_bit_width) + lvl2_match for lvl2_match in
                                        list(self._level_2[lvl1_match].match(inf_value, inf_ignore))],
                    lvl1_matches)
            )
        else:
            return lvl1_matches

    def bit_length(self):
        """
        Number of bits needed to encode the stored integer sequence.
        """
        return self._level_1.bit_length() + sum(ef2.bit_length() for ef2 in (self._level_2.values()))

    def __len__(self) -> int:
        """
        Return number of integers stored in the Elias-Fano structure.
        """
        return self._n

    def __getitem__(self, k: int) -> int:
        """
        Return k-th integer stored in this Elias-Fano structure using subscript operator.
        """
        return self.select(k)

    def __iter__(self) -> Iterator[int]:
        """
        Support for __iter__ and next
        """
        if bool(self._level_2):
            # recursively prepend lvl-1 element (aka sup) to each of its respective lvl-2 elements (aka inf)

            return chain.from_iterable(
                map(lambda p: [(p[0] << (int(math.log2(self._u)) - self._b)) + v for v in list(iter(p[1]))],
                    iter(self._level_2.items())))
        else:
            return iter(self._level_1)

    def to_bytes(self):
        """

        """
        _level_1_byte_count, _level_1_bytes = self._level_1.to_bytes()
        _level_2_count = len(self._level_2)

        rep_size_n = max(1, math.ceil(self._n.bit_length() / 8.0))
        rep_size_level_1_byte_count = max(1, math.ceil(_level_1_byte_count.bit_length() / 8.0))
        rep_size_level_2_count = max(1, math.ceil(_level_2_count.bit_length() / 8.0))

        if _level_2_count:
            _level_2_byte_counts, _level_2_bytes = unzip(map(lambda lvl2: lvl2.to_bytes(),
                                                             map(itemgetter(1),
                                                                 sorted(self._level_2.items()))
                                                             )
                                                         )
            rep_size_max_level_2_byte_count = max(1, math.ceil(max(_level_2_byte_counts).bit_length() / 8.0))
        else:
            _level_2_byte_counts, _level_2_bytes, rep_size_max_level_2_byte_count = 0

        ba = []

        # HEADER
        #  Index type:                     1 byte                                      EliasFano: 0, MultiLevelEliasFano: 1
        #  rep_size_n:                     1 byte                                      #bytes for _n
        #  rep_size_level_1_byte_count:    1 byte                                      #bytes for _level_1_byte_count
        #  rep_size_level_2_count:         1 byte                                      #bytes to store _level_2_count
        #  rep_size_max_level_2_byte_count 1 byte                                      #bytes for the largest _level_2_byte_count

        #  _n:                             (rep_size_n) bytes                          #elements in sequence
        #  _max_bits:                      1 byte                                      allows reconstruction of _u
        #  _d:                             1 byte                                      depth
        #  _b:                             1 byte                                      #bits handled at lvl_1

        #  _level_1_byte_count:            (rep_size_level_1_byte_count) bytes         #bytes to store lvl1 index

        #  _level_2_count:                 (rep_size_level_2_count) bytes              #lvl2_indices
        #  _level_2_byte_count:            (rep_size_max_level_2_byte_count) bytes     #bytes for each lvl2 index
        #  ...
        #  _level_2_byte_count:            (rep_size_max_level_2_byte_count) bytes     #bytes for each lvl2 index

        ba.append(int(1).to_bytes(1, 'little', signed=False))
        ba.append(rep_size_n.to_bytes(1, 'little', signed=False))
        ba.append(rep_size_level_1_byte_count.to_bytes(1, 'little', signed=False))
        ba.append(rep_size_level_2_count.to_bytes(1, 'little', signed=False))
        ba.append(rep_size_max_level_2_byte_count.to_bytes(1, 'little', signed=False))

        ba.append(self._n.to_bytes(rep_size_n, 'little', signed=False))
        ba.append(int(math.log2(self._u)).to_bytes(1, "little", signed=False))
        ba.append(self._d.to_bytes(1, 'little', signed=False))
        ba.append(self._b.to_bytes(1, 'little', signed=False))

        ba.append(_level_1_byte_count.to_bytes(rep_size_level_1_byte_count, 'little', signed=False))
        ba.append(_level_2_count.to_bytes(rep_size_level_2_count, 'little', signed=False))

        if _level_2_count:
            for _level_2_byte_count in _level_2_byte_counts:
                ba.append(_level_2_byte_count.to_bytes(rep_size_max_level_2_byte_count, 'little', signed=False))

        # DATA
        #   _level_1_bytes:         (_level_1_byte_counts) bytes
        #   _level_2_bytes:         sum(_level_2_byte_count) bytes
        ba.append(_level_1_bytes)

        if _level_2_count:
            for _level_2_byte in _level_2_bytes:
                ba.append(_level_2_byte)

        return len(b''.join(ba)), b''.join(ba)

    @staticmethod
    def from_bytes(index_bytes: bytearray):

        bytes_iter = iter(index_bytes)

        # HEADER
        #  Index type:                     1 byte                                      EliasFano: 0, MultiLevelEliasFano: 1
        #  rep_size_n:                     1 byte                                      #bytes for _n
        #  rep_size_level_1_byte_count:    1 byte                                      #bytes for _level_1_byte_count
        #  rep_size_level_2_count:         1 byte                                      #bytes to store _level_2_count
        #  rep_size_max_level_2_byte_count 1 byte                                      #bytes for the largest _level_2_byte_count

        #  _n:                             (rep_size_n) bytes                          #elements in sequence
        #  _max_bits:                      1 byte                                      allows reconstruction of _u
        #  _d:                             1 byte                                      depth
        #  _b:                             1 byte                                      #bits handled at lvl_1

        #  _level_1_byte_count:            (rep_size_level_1_byte_count) bytes         #bytes to store lvl1 index
        #  _level_2_count:                 (rep_size_level_2_count) bytes              #lvl2_indices

        #  _level_2_byte_count:            (rep_size_max_level_2_byte_count) bytes     #bytes for each lvl2 index
        #  ...
        #  _level_2_byte_count:            (rep_size_max_level_2_byte_count) bytes     #bytes for each lvl2 index

        _index_type = int.from_bytes(take(1, bytes_iter), 'little', signed=False)
        rep_size_n = int.from_bytes(take(1, bytes_iter), 'little', signed=False)
        rep_size_level_1_byte_count = int.from_bytes(take(1, bytes_iter), 'little', signed=False)
        rep_size_level_2_count = int.from_bytes(take(1, bytes_iter), 'little', signed=False)
        rep_size_max_level_2_byte_count = int.from_bytes(take(1, bytes_iter), 'little', signed=False)

        _n = int.from_bytes(take(rep_size_n, bytes_iter), 'little', signed=False)
        _max_bits = int.from_bytes(take(1, bytes_iter), 'little', signed=False)
        _d = int.from_bytes(take(1, bytes_iter), 'little', signed=False)
        _b = int.from_bytes(take(1, bytes_iter), 'little', signed=False)

        _level_1_byte_count = int.from_bytes(take(rep_size_level_1_byte_count, bytes_iter), 'little', signed=False)
        _level_2_count = int.from_bytes(take(rep_size_level_2_count, bytes_iter), 'little', signed=False)

        if _level_2_count:
            _level_2_byte_counts = [int.from_bytes(take(rep_size_max_level_2_byte_count, bytes_iter), 'little', signed=False) for _ in range(_level_2_count)]

        # DATA
        #   _level_1_bytes:         (_level_1_byte_counts) bytes
        #   _level_2_bytes:         sum(_level_2_byte_count) bytes
        _level_1_bytes = take(_level_1_byte_count, bytes_iter)

        if _level_2_count:
            _level_2_bytes = [take(cnt, bytes_iter) for cnt in _level_2_byte_counts]

        # TODO: implement appropriate constructor
        mef = MultiLevelEliasFano([0], 1)
        mef._n = _n
        mef._b = _b
        mef._d = _d
        mef._u = 2 ** _max_bits
        mef._level_1 = EliasFano.from_bytes(_level_1_bytes)
        mef._level_2 = {}

        if _level_2_count:
            _level_2_indices = [MultiLevelEliasFano.from_bytes(x) if x[0] == 1 else EliasFano.from_bytes(x)
                                for x in _level_2_bytes]

            for label, index in zip(iter(mef._level_1), _level_2_indices):
                mef._level_2[label] = index

        else:
            mef._level_2 = {}

        return mef

    def to_file(self, file_path: str):
        """
        Save given Elias-Fano structure to disk.
        :param file_path: file to hold the EF structure
        :return: None
        """
        byte_size, all_bytes = self.to_bytes()

        with open(file_path, 'wb') as file:
            file.write(all_bytes)

        logging.info("Serialized multi-level EF index using %s bytes." % byte_size)

    @staticmethod
    def from_file(file_path: str):
        """
        Load multi-level Elias-Fano structure from disk.
        :param file_path: file storing an Elias-Fano structure
        :return: an Elias-Fano structure
        """
        import os

        assert os.path.exists(file_path) and os.path.isfile(file_path), IOError("File path invalid or does not exist.")

        file = open(file_path, "rb")
        index_bytes = bytearray(file.read(os.path.getsize(file_path)))
        file.close()

        return MultiLevelEliasFano.from_bytes(index_bytes)
