# kbbq

[![codecov](https://codecov.io/gh/adamjorr/kbbq/branch/master/graph/badge.svg)](https://codecov.io/gh/adamjorr/kbbq)
[![CircleCI](https://circleci.com/gh/adamjorr/kbbq.svg?style=svg)](https://circleci.com/gh/adamjorr/kbbq)
[![Documentation Status](https://readthedocs.org/projects/kbbq/badge/?version=latest)](https://kbbq.readthedocs.io/en/latest/?badge=latest)


k-mer based base quality recalibration

## a helpful k-mer diagram

```
               +
               |
               v
seqlen=7    0123456
nkmers=5    ATGCATG
ksize=3   +0XXX
          +1 XXX
           2  XXX
          +3   XXX
          +4    XXX
```
