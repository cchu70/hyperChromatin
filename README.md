# hyperChromatin

APPMTH 220 course project - application of hyperbolic spaces for representing single cell lineage using chromatin accessibility.

## Interactive plot demo
See notebook 10

<img src="https://github.com/cchu70/hyperChromatin/raw/main/interactive_poincare_ball.gif" width="400"/>

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

# umap issues
pip install --upgrade tbb
pip install --upgrade --force-reinstall scanpy
pip install tbb-devel
# pip install --upgrade --force-reinstall numba
pip install numba==0.56.2
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
