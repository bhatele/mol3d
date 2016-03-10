mol3d
=====

Authors:
```
  Abhinav Bhatele <bhatele@llnl.gov>
  Edgar Solomonik
```

### Build

1. Run the copy_namd.sh script (specify a correct path for namd inside the
script) in the src_namd folder
NOTE: NAMD source from 2010-04-01 will definitely work
2. Go to the compile folder and type:
```
make
```

### Run

To run specify three parameters:
1. A NAMD configuration file (such as apoa1.namd)
2. The branching factor (default set to 2), and
3. A load balancing strategy

Example:
```bash
./charmrun +p8 ./mol3d apoa1.namd 8 +balancer DummyLB
```
Run from the same folder as namd parameter files (need a .js format file)
