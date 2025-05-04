# hyperChromatin

APPMTH 220 course project - application of hyperbolic spaces for representing single cell lineage using chromatin accessibility.


# Installation

## SIMBA
```bash
# use pyenv
python3 -m venv env/hyperChrom_pyenv
source env/hyperChrom_pyenv/bin/activate
pip install git+https://github.com/huidongchen/simba
pip install ipykernel ipython
python -m ipykernel install --name=hyperChrom_pyenv
# restart kernel
source env/hyperChrom_pyenv/bin/activate
pip install --upgrade attrs

# simba pbg
cd src
git submodule addgit@github.com:pinellolab/simba_pbg.git
cd simba_pbg
pip install -e .

pip install scanpy
```

## Poincare map
```bash
cd src
git submodule add git@github.com:facebookresearch/PoincareMaps.git

source env/hyperChrom_pyenv/bin/activate
pip install fastdtw
```

## scDHMap

```bash
cd src
git submodule add git@github.com:ttgump/scDHMap.git
```

## scPhere

```bash
cd src
git submodule add git@github.com:klarman-cell-observatory/scPhere.git
```
