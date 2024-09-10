
# This script installs Miniconda and creates a new conda environment called pyquest
#!/bin/bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh

chmod +x ~/miniconda.sh

bash ~/miniconda.sh -b -p $HOME/miniconda

echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc

source $HOME/.bashrc

conda init bash

conda create -n pyquest python=3.9
conda activate pyquest
git clone -b develop --recursive https://github.com/rrmeister/pyQuEST
pip install ./pyQuEST
pip install -U numpy==1.26.4 matplotlib pandas tqdm
