# pyEliasFano
[![DOI](https://zenodo.org/badge/367291041.svg)](https://zenodo.org/badge/latestdoi/367291041)

pyEliasFano offers **quasi-succinct** representations for monotone non-decreasing sequences of n integers from 
a universe [0 . . . m). 

We currently support the following variants of Elias-Fano indexing:
* ``pyEliasFano.EliasFano``: the classical Elias-Fano representation occupying occupying 2*n+n*ceil(log2(m)/n) bits 
* ``pyEliasFano.UniformlyPartitionedEliasFano``: an uniformly-partitioned Elias-Fano representation 

All variants support the following operations:
- ``select(i)``: fast access to the ``i``-th element of the integer sequence,
- ``rank(x)``: queries the index position of ``x`` iff stored within the given Elias-Fano structure 
- ``nextGEQ(x)``: fast access to the smallest integer of the sequence that is greater or equal than ``x``
- ``nextLEQ(x)``: fast access to the largest integer of the sequence that is smaller or equal than ``x``

## Installation
Install from PyPi using
```bash
pip install pyEliasFano
```

## Usage
```python
from pyEliasFano import EliasFano, UniformlyPartitionedEliasFano
```
imports the module.

```python
integers = sorted([123, 1343, 2141, 35312, 4343434])
ef0 = EliasFano(integers)
ef1 = UniformlyPartitionedEliasFano(integers)
```
creates a classical Elias-Fano structure ``ef0`` as well as an uniformly-partitioned Elias-Fano structure ``ef1`` for the **sorted** ``integers`` sequence. 

### Access
The ``i``th element from the original ``integers`` sequence can be retrieved from an Elias-Fano structure ``ef`` using its ``select(i)`` method
```python
ef0.select(3)
>>> 35312

ef0.select(0)
>>> 123
```
or using subscript operator
```python
ef1[3]
>>> 35312
```

An Elias-Fano structure ``ef`` is also iterable. 

You can easily loop through the stored elements stored 
```python
ef_iter = iter(ef0)

next(ef_iter)
>>> 123

next(ef_iter)
>>> 1343   
```
or return all stored elements at once
```python
list(iter(ef0))
>>> [123, 1343, 2141, 35312, 4343434]
```
As a side note, the following assertion will always hold:
```python
assert [ef.select(ef.rank(v)) for v in integers] == integers
```

### Rank
Given an integer ``x``, we can query the index position of ``x`` within an Elias-Fano structure ``ef`` using its ``rank(x)`` method.

For example,
```python
ef0.rank(4343434)
>>> 4

ef1.rank(123)
>>> 0
```

As a side note, the following assertion will always hold:
```python
assert [ef.rank(ef.select(i)) for i in range(len(integers))]
```

### nextGEQ
Given an integer ``x``, we can query the smallest integer stored within an Elias-Fano structure ``ef`` that is larger than or equal to ``x`` using the ``nextGEQ(x)``method.

For example,
```python
ef0.nextGEQ(1345)
>>> 2141

ef0.nextGEQ(4343434)
>>> 4343434

ef1.nextGEQ(2)
>>> 123
```

### nextLEQ
Given an integer ``x``, we can query the largest integer stored within an Elias-Fano structure ``ef`` that is smaller than or equal to ``x`` using the ``nextLEQ(x)``method.

For example,
```python
ef0.nextLEQ(4343420)
>>> 35312

ef0.nextLEQ(123)
>>> 123
```

# Citation
```bibtex
@misc{rmrschub_2021_pyEliasFano,
    author       = {René Schubotz},
    title        = {{pyEliasFano: Quasi-succinct represenations for monotone non-decreasing sequences of integers.}},
    month        = may,
    year         = 2021,
    doi          = {10.5281/zenodo.4774741},
    version      = {0.0.6},
    publisher    = {Zenodo},
    url          = {https://github.com/rmrschub/pyEliasFano}
    }
```

# License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
