#!/usr/bin/env python3

import h5py
import numpy as np

if __name__ == "__main__":
    a = np.random.random(size=(16,16,16)).astype(np.float32)
    f = h5py.File("inputs_float.h5", "w")
    grp = f.create_group("level_0")
    dset = grp.create_dataset("dark_matter_density", data=a)
    f.flush()