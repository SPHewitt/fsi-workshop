# Install instructions for AWS Ubuntu 16.04

## Foam-extend 4.0
Follow the below instructions to install Foam-extend 4.0.

```bash
sudo apt-get update

sudo apt-get install git-core build-essential binutils-dev cmake flex \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz python python-dev

cd ~
mkdir foam
cd foam
git clone git://git.code.sf.net/p/foam-extend/foam-extend-4.0 foam-extend-4.0
cd ~/foam/foam-extend-4.0
source etc/bashrc
echo "alias fe40='source \$HOME/foam/foam-extend-4.0/etc/bashrc'" >> $HOME/.bashrc
export QT_BIN_DIR=/usr/bin/
echo "export QT_BIN_DIR=$QT_BIN_DIR" >> etc/prefs.sh
./Allwmake.firstInstall -j 4
```

## ParaFEM

```bash
cd ~
git clone https://github.com/leemargetts/ParaFEM.git
cd ParaFEM/parafem
export PATH=$WM_THIRD_PARTY_DIR/packages/openmpi-1.8.8/platforms/linux64GccDPOpt/bin:$PATH
sed -i "s:/usr/bin/mpif90:mpif90:" build/linuxdesktop.inc
./make-parafem MACHINE=linuxdesktop
```
Errors: No underlying compiler was specified in the wrapper compiler data file
(e.g., mpicc-wrapper-data.txt)

Solution: Install a version of gFortran on the Machine and reinstall the OpenMPI package in OpenFOAM
sudo apt-get install gfortran

## Fluid-structure interaction library
```bash
mkdir -p $FOAM_RUN
cd $FOAM_RUN/..
```


