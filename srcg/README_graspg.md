# GRASPG:  An extension of the GRASP-2018 package, using Configuration State Function Generators

GRASPG is an extension of the Grasp2018 package [Comput. Phys. Commun. 237 (2019) 184187]
based on Configuration State Function Generators (CSFGs). The ideas could be seen in
[Comput. Phys. Commun. 283 (2023) 108562] and [Comput. Phys. Commun. XXX (20XX) XXXXXX],
in the latter the long write-up can be found.

## Installation

1. Go to the main directory grasp-master of the installed Grasp2018 package. 
    Run 
              source make-environment_xxx
    , where xxx is the compiler name. The Grasp2018 environment variables are now set.
2. Copy the graspg package, downloaded from GitHub or the CPC website, to the grasp-master directory. 
    Untar the package, the directories srcg and graspgtest will now appear.
3. In the grasp-master/srcg directory, execute the installation by issuing the commands
             make clean & make
   This will generate static library archives that reside in the grasp-master/lib directory. 
   It will also generate executable program files that reside in the grasp-master/bin directory. 
   Graspg is now fully merged with the Grasp2018 package.

## About GRASPG

This version of GRASPG is an extension of the Grasp2018 (https://github.com/compas/grasp) package.
In this version, the original rci_mpi program is modified by using Configuration State Function Generators. 
Then this package should be referred to together with the original GRASP package, 
i.e. [Comput. Phys. Commun. 237 (2019) 184187](https://doi.org/10.1016/j.cpc.2018.10.032). 

The ideas of this version of GRASPG have been published in:
> Yan TingLi, Kai Wang, Ran Si, Michel Godefroid, Gediminas Gaigalas, ChongYang Chen, P. Jönsson, 
> "Reducing the computational load - atomic multiconfiguration calculations based on configuration state function generators",
> Computer Physics Communications, 283, 108562 (2023),
> https://doi.org/10.1016/j.cpc.2022.108562 

The long write-up and codes have been published in :
> Ran Si, Yanting Li, Kai Wang, Chongyang Chen, Gediminas Gaigalasc, Michel Godefroidd, Per Jönsson，
> "Graspg -- An extension to Grasp2018 based on Configuration State Function Generators"，
> Computer Physics Communications, XXX, XXXXXX (20XX)

Development of this package was performed largely by:
|                                          | email                         |
| ------------------------- | ------------------------------|
| Chongyang Chen           | chychen@fudan.edu.cn          |
| Gediminas Gaigalas       | Gediminas.Gaigalas@tfai.vu.lt |
| Per Jönsson                    | per.jonsson@mau.se            |

Please contact the repository manager should you have any questions with regards
to bugs or the general development procedure. Contact the leading developer for 
specific questions related to a certain code.

## Structure of the Package

The package has the structure shown below where executables, after successful
compilation, reside in the `bin` directory. Compiled libraries are in the `lib`
directory. Scripts for example runs and case studies are in folders under
`graspgtest`. Source code is in the `srcg` directory and divided into applications
in the `appl` directory, libraries in the `lib` directory and tools in the
`tool` directory.

```
  |-graspgtest
  |-----script1
  |-----script2
  |-----script3
  |-----script3_ZF
  |-----script4
  |-----script5
  |-lib
  |---lib9290_csfg
  |---libdvd90mpi
  |---libme_csfg
  |---libmod_csfg
  |---librang90_csfg
  |-src
  |---appl
  |-----jj2lsj90_csfg
  |-----rangular90_csfg_mpi
  |-----rci90_block_csfg_mpi
  |-----rci90_csfg_mpi
  |-----rcsfggenerate90_csfg
  |-----rcsfginteract90_csfg
  |-----rcsfgzerofirst90_csfg
  |-----rmcdhf90_csfg_mpi  
  |-----rmixaccumulate_csfg
  |---tool
