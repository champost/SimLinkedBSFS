## `SimLinkedBSFS`

### Description
`SimLinkedBSFS` is a C++ wrapper for the popular command line software for simulating population genetic samples [`ms`](http://home.uchicago.edu/rhudson1/source/mksamples.html). The role of `SimLinkedBSFS` is to generate a `bSFS` (or blockwise Site Frequency Spectrum) from the `ms` output which is already in blocks but in a sequence format. `SimLinkedBSFS` is most useful when simulating large stretches of sequence (say 1Mbp in size) but when the `bSFS` needs to be generated from much smaller blocks (say 2Kb in size).

### Install instructions
It is easiest to build a `SimLinkedBSFS` binary under all flavours of Linux. `SimLinkedBSFS` requires the GNU Compiler Collection ([`gcc`](https://gcc.gnu.org/)). If you do not have `gcc`, you can use your OS specific package handling utility (under Debian/Ubuntu) :

    sudo apt-get install build-essential

1. Once `gcc` is available, run

       wget https://github.com/champost/SimLinkedBSFS/archive/master.tar.gz

2. Untar the archive using 

       tar -xzf master.tar.gz

3. Change directory

        cd SimLinkedBSFS-master

4. Make the compile script executable

        chmod +x compile.sh

5. Compile `SimLinkedBSFS`

        ./compile.sh

If you want `SimLinkedBSFS` to be accessible from anywhere, such as your data folder, you might want to

    cp SimLinkedBSFS ~/bin

This ensures that you can execute the program by specifying `SimLinkedBSFS ...` instead of `./SimLinkedBSFS ...` from the installation folder. This holds only if `~/bin` exists and is part of your `$PATH` environment variable.

### Usage
The commands for `SimLinkedBSFS` are almost entirely in the [`ms`](http://home.uchicago.edu/rhudson1/source/mksamples.html) format. For example we would like to simulate a 2 population scenario consisting of 5 and 2 haploid samples respectively. The command below will simulate 100000 sequences of length 4000 bp.

    ms 7 100000 -t 1 -r 0.5 4001 -I 2 5 2 -m 1 2 1 -ej 1.25 1 2

Now let us assume that we want to create a `bSFS` consisting of 1000 bp blocks. We start by creating a file named `popconfig.txt` containing the number of samples per population. Here it would be :

    echo 5 2 > popconfig.txt

Next, we can run `SimLinkedBSFS` as follows :

    SimLinkedBSFS 7 100000 -t 1 -r 0.5 4001 -I 2 5 2 -m 1 2 1 -ej 1.25 1 2 A B C D > my_two_pop_bSFS.txt

Here `A, B, C, D` are what differentiate it from the `ms` command line.

    A : a value z (0 < z <= 1) giving the block size (here it would be 4000 X z = 1000, so z = 0.25)
    B : number of mutations per branch class or Kmax (please refer to the bSFS literature). 0 means no Kmax.
    C : number of simulated populations (here it will be 2)
    D : a lowercase "f" or "u", saying you want the folded or unfolded bSFS respectively

Finally, the command line for generating a folded bSFS from 1000 bp blocks and no Kmax would be :

    SimLinkedBSFS 7 100000 -t 1 -r 0.5 4001 -I 2 5 2 -m 1 2 1 -ej 1.25 1 2 0.25 0 2 f > my_two_pop_bSFS.txt









