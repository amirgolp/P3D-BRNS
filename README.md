# About P3D-BRNS

P3D-BRNS is an open-source modelling package designed for simulating Multiphae, Multi-component reactive transport in &mu;-CT image of any porous strcuture. 
The "lib" folder contains the code for tracking interface proeprties (Raeini et al 2016), while all others contain the relevant solver and 
test cases used in the main manuscript. 

# Instructions for running P3D-BRNS locally:

The source code is developed in a linux-based environemnt (Linux Ubuntu 18.04) using OpenFOAM-4.0 extended version.

## Step 1: 
First step is to install a working version of the OpenFOAM-4.0x. For this open a shell and type the followings:

1.
```sh 
sudo apt-get update
```
2.
```sh 
sudo apt-get install git-core build-essential binutils-dev cmake flex zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libxt-dev rpm mercurial graphviz python python-dev gcc-5 g++-5
```
4.
```sh 
cd ~
```
5.
```sh 
mkdir foam
```
6.
```sh 
cd foam
```
7.
```sh 
git clone git://git.code.sf.net/p/foam-extend/foam-extend-4.0 foam-extend-4.0
```
8.
```sh 
cd ~/foam/foam-extend-4.0
```
9.
```sh 
echo export WM_THIRD_PARTY_USE_BISON_27=1 >> etc/prefs.sh
```
10.
```sh 
echo export QT_SELECT=qt4 >> etc/prefs.sh
```
11.
```sh 
echo "export WM_CC='gcc-5'" >> etc/prefs.sh
```
12.
```sh 
echo "export WM_CXX='g++-5'" >> etc/prefs.sh
```
13.
```sh 
echo "export QT_BIN_DIR=/usr/bin/" >> etc/prefs.sh   
```
14.
```sh 
source etc/bashrc
```
15.
```sh 
echo "alias fe40='source \$HOME/foam/foam-extend-4.0/etc/bashrc'" >> $HOME/.bashrc
```
16.
```sh 
sed -i -e 's=rpmbuild --define=rpmbuild --define "_build_id_links none" --define=' ThirdParty/tools/makeThirdPartyFunctionsForRPM
```
17.
```sh 
sed -i -e 's/gcc/\$(WM_CC)/' wmake/rules/linux64Gcc/c
```
18.
```sh 
sed -i -e 's/g++/\$(WM_CXX)/' wmake/rules/linux64Gcc/c++ 
```
19.
```sh 
./Allwmake.firstInstall
```


You can test for successfull installation by typing the following in the shell:
```sh
  interFoam -help
```
which should produce no error.

## Step 2
Next, the library @ "lib/interfacePropertiesAM" needs to be compiled. For this, navigate to this folder and in the shell, type:
``` sh
wmake libso
```

## Step 3:
To run each test case, first the corresponding solver needs to be compiled. Navigate to solver folder and in the shell, type:
``` sh
wmake
```
Now navigate to the relative test-case fodler and in the shell, type the name of the solver.




