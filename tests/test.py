import os
from pyEliasFano import EliasFano, save, load

if __name__ == "__main__":

    # create an Elias-Fano structure
    values = sorted([123, 1343, 2141, 35312, 4343434])
    ef = EliasFano(values)

    # store EF structure to disk
    save(ef, "example_elias_fano.pickle")

    # select(k)
    assert ef.select(0) == 123
    assert ef.select(1) == 1343
    assert ef.select(2) == 2141
    assert ef.select(3) == 35312
    assert ef.select(4) == 4343434

    # rank(x)
    assert ef.rank(123) == 0
    assert ef.rank(1343) == 1
    assert ef.rank(2141) == 2
    assert ef.rank(35312) == 3
    assert ef.rank(4343434) == 4

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
    assert ef.nextLEQ(1343) == 123
    assert ef.nextLEQ(2141) == 1343
    assert ef.nextLEQ(35312) == 2141
    assert ef.nextLEQ(4343434) == 4343434
    assert ef.nextLEQ(500000000000) == 4343434

    del ef

    # load EF structure from disk
    ef2 = load("example_elias_fano.pickle")

    # select(k)
    assert ef2.select(0) == 123
    assert ef2.select(1) == 1343
    assert ef2.select(2) == 2141
    assert ef2.select(3) == 35312
    assert ef2.select(4) == 4343434

    # rank(x)
    assert ef2.rank(123) == 0
    assert ef2.rank(1343) == 1
    assert ef2.rank(2141) == 2
    assert ef2.rank(35312) == 3
    assert ef2.rank(4343434) == 4

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
    assert ef2.nextLEQ(1343) == 123
    assert ef2.nextLEQ(2141) == 1343
    assert ef2.nextLEQ(35312) == 2141
    assert ef2.nextLEQ(4343434) == 4343434
    assert ef2.nextLEQ(500000000000) == 4343434

    # cleanup
    try:
        os.remove("example_elias_fano.pickle")
    except OSError:
        pass