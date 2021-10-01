# pyEliasFano
[![DOI](https://zenodo.org/badge/367291041.svg)](https://zenodo.org/badge/latestdoi/367291041)

pyEliasFano offers a **quasi-succinct** represenation for a monotone non-decreasing sequence of n integers from 
the universe [0 . . . m) occupying 2*n+n*ceil(log2(m/n)) bits.

![](https://www.antoniomallia.it/uploads/Elias-Fano.png)
*by Antonio Mallia, 2018*

It supports the following operations:
- **select(k)**: nearly constant time access to the k-th element,
- **rank(x)**: access to the index within the structure for given integer x.
- **nextGEQ(x)**: fast access to the smallest integer of the sequence that is greater or equal than a given x
- **nextLEQ(x)**: fast access to the largest integer of the sequence that is smaller or equal than a given x

## Installation
```bash
pip install pyEliasFano
```

## Usage
```python
from pyEliasFano import EliasFano
```
imports the module.

```python
integers = sorted([123, 1343, 2141, 35312, 4343434])
ef = EliasFano(integers)
```
creates an Elias-Fano structure for the sorted ``integers`` list. 

### Access
A single element stored within the EliasFano structure can be accessed using the ``select`` method
```python
ef.select(3)
```
or the subscript operator
```python
ef[3]
```
{% hint style="info" %}
Indexing in an EliasFano structure is one-based!
% endhint %}

The EliasFano structure is also iterable. You can easily loop through the elements stored in an EliasFano structure
```python
for v in iter(ef):
    print(v)    
```
or return all elements at once. 
```python
list(iter(ef))
```

### Rank
```python
ef.rank(4343434)
```
returns the index position for the given integer iff stored within the Elias-Fano structure. 
Here, we get ``5``.

### nextGEQ
```python
ef.nextGEQ(1345)
```
returns the smallest integer stored in the Elias-Fano structure that is larger than or equal to the given integer. 
### nextLEQ
```python
ef.nextLEQ(4343420)
```
returns the largest integer stored in the Elias-Fano structure that is smaller than or equal to the given integer.
Here, we get ``35312``. 

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
