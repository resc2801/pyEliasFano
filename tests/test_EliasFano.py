import os
from pyEliasFano import EliasFano, save, load

if __name__ == "__main__":

    # create an Elias-Fano structure
    values = sorted([123, 1343, 2141, 35312, 4343434])
    ef = EliasFano(values)

    # store EF structure to disk
    save(ef, "example_elias_fano.pickle")

    assert [ef.select(ef.rank(v)) for v in values] == values
    assert [ef.rank(ef.select(i)) for i in range(len(values))] == list(range(len(values)))

    # nextGEQ(x)
    assert ef.nextGEQ(0) == 123
    assert ef.nextGEQ(75) == 123
    assert ef.nextGEQ(123) == 123
    assert ef.nextGEQ(231) == 1343
    assert ef.nextGEQ(1343) == 1343
    assert ef.nextGEQ(1750) == 2141
    assert ef.nextGEQ(2141) == 2141
    assert ef.nextGEQ(35312) == 35312
    assert ef.nextGEQ(39312) == 4343434
    assert ef.nextGEQ(4343434) == 4343434

    # nextLEQ(x)
    assert ef.nextLEQ(123) == 123
    assert ef.nextLEQ(1343) == 1343
    assert ef.nextLEQ(2141) == 2141
    assert ef.nextLEQ(35312) == 35312
    assert ef.nextLEQ(4343434) == 4343434
    assert ef.nextLEQ(500000000000) == 4343434

    del ef

    # load EF structure from disk
    ef2 = load("example_elias_fano.pickle")

    assert [ef2.select(ef2.rank(v)) for v in values] == values
    assert [ef2.rank(ef2.select(i)) for i in range(1, len(values)+1)] == list(range(1, len(values)+1))

    # nextGEQ(x)
    assert ef2.nextGEQ(0) == 123
    assert ef2.nextGEQ(75) == 123
    assert ef2.nextGEQ(123) == 123
    assert ef2.nextGEQ(231) == 1343
    assert ef2.nextGEQ(1343) == 1343
    assert ef2.nextGEQ(1750) == 2141
    assert ef2.nextGEQ(2141) == 2141
    assert ef2.nextGEQ(35312) == 35312
    assert ef2.nextGEQ(39312) == 4343434
    assert ef2.nextGEQ(4343434) == 4343434

    # nextLEQ(x)
    assert ef2.nextLEQ(123) == 123
    assert ef2.nextLEQ(191) == 123
    assert ef2.nextLEQ(1343) == 1343
    assert ef2.nextLEQ(1743) == 1343
    assert ef2.nextLEQ(2141) == 2141
    assert ef2.nextLEQ(8141) == 2141
    assert ef2.nextLEQ(35312) == 35312
    assert ef2.nextLEQ(353120) == 35312
    assert ef2.nextLEQ(4343434) == 4343434
    assert ef2.nextLEQ(500000000000) == 4343434

    # cleanup
    try:
        os.remove("example_elias_fano.pickle")
    except OSError:
        pass