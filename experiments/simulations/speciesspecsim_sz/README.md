# Compilation Instructions

## Prerequisites

Compilation requires GNU Make and the GNU Scientific Library (GSL).

Ensure that GSL is installed and its shared libraries are accessible.

## Compilation

First, set the `LD_LIBRARY_PATH` environment variable to include the GSL library directory:

```bash
export LD_LIBRARY_PATH=/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/gsl/lib:$LD_LIBRARY_PATH
```

Replace the path above with the location of your GSL installation if necessary.

Then compile the program:

```bash
make
```

**Note:** The `LD_LIBRARY_PATH` setting may also be required at runtime to locate the GSL shared libraries.