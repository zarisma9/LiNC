#!/bin/bash

set -euo pipefail

# ========================================================
# User-defined settings
# ========================================================
# Check the restraint mask for your system.
# Modify this mask if additional solvent/ion/residue names are present.
RESTRAINT_MASK="!:WAT,Na+,K+,Cl-"

echo ">>> Generating AMBER input files (Steps 1-7)..."

# ========================================================
# 1. Minimization 1
#    Solvent minimization with strong restraints on solute
# ========================================================

cat <<EOF > 01_Min1.in
Minimize solvent; restrain solute strongly
 &cntrl
  imin=1, maxcyc=3000, ncyc=1000,
  ntb=1, cut=10.0,
  ntf=2, ntc=2,
  ntpr=100, ntwr=1000,
  ntr=1, restraint_wt=10.0,
  restraintmask='${RESTRAINT_MASK}',
 /
EOF

# ========================================================
# 2. Minimization 2
#    Whole-system minimization with weak restraints
# ========================================================

cat <<EOF > 02_Min2.in
Whole system minimization with weak restraints
 &cntrl
  imin=1, maxcyc=10000, ncyc=5000,
  ntb=1, cut=10.0,
  ntf=2, ntc=2,
  ntpr=100, ntwr=2000,
  ntr=1, restraint_wt=2.0,
  restraintmask='${RESTRAINT_MASK}',
 /
EOF

# ========================================================
# 3. Heating
#    NVT heating from 0 K to 310 K
# ========================================================

cat <<EOF > 03_Heat_NVT.in
Heat from 0 K to 310 K at NVT with restraints on solute
 &cntrl
  imin=0, ntx=1, irest=0,
  ntc=2, ntf=2, tol=1.0e-7,
  dt=0.001, nstlim=25000,
  ntt=3, gamma_ln=1.0, ig=-1,
  ntpr=100, ntwr=10000, ntwx=100,
  nmropt=1,
  ntb=1, ntp=0, cut=10.0, ioutfm=1, ntxo=2,
  ntr=1, restraint_wt=5.0,
  restraintmask='${RESTRAINT_MASK}',
 /
 &wt type='TEMP0', istep1=0, istep2=25000, value1=0.0, value2=310.0 /
 &wt type='END' /
EOF

# ========================================================
# 4. Equilibration 1
#    NPT equilibration with restraint weight 2.0
# ========================================================

cat <<EOF > 04_Equil_NPT_1.in
Equilibration NPT with restraints, approximately 200 ps
 &cntrl
  imin=0, ntx=5, irest=1,
  ntc=2, ntf=2, tol=1.0e-7,
  dt=0.002, nstlim=100000,
  ntt=3, gamma_ln=1.0, temp0=310.0, ig=-1,
  ntpr=1000, ntwr=50000, ntwx=1000,
  ntb=2, cut=10.0, ioutfm=1, ntxo=2,
  ntp=1, taup=1.0, barostat=2,
  ntr=1, restraint_wt=2.0,
  restraintmask='${RESTRAINT_MASK}',
 /
EOF

# ========================================================
# 5. Equilibration 2
#    NPT equilibration with weaker restraint weight 0.5
# ========================================================

cat <<EOF > 05_Equil_NPT_2.in
Equilibration NPT with weaker restraints, approximately 200 ps
 &cntrl
  imin=0, ntx=5, irest=1,
  ntc=2, ntf=2, tol=1.0e-7,
  dt=0.002, nstlim=100000,
  ntt=3, gamma_ln=1.0, temp0=310.0, ig=-1,
  ntpr=1000, ntwr=50000, ntwx=1000,
  ntb=2, cut=10.0, ioutfm=1, ntxo=2,
  ntp=1, taup=1.0, barostat=2,
  ntr=1, restraint_wt=0.5,
  restraintmask='${RESTRAINT_MASK}',
 /
EOF

# ========================================================
# 6. Equilibration 3
#    Unrestrained NPT equilibration
# ========================================================

cat <<EOF > 06_Equil_NPT_3.in
Equilibration NPT without restraints, approximately 500 ps
 &cntrl
  imin=0, ntx=5, irest=1,
  ntc=2, ntf=2, tol=1.0e-7,
  dt=0.002, nstlim=250000,
  ntt=3, gamma_ln=1.0, temp0=310.0, ig=-1,
  ntpr=2500, ntwr=125000, ntwx=2500,
  ntb=2, cut=10.0, ioutfm=1, ntxo=2,
  ntp=1, taup=1.0, barostat=2,
  ntr=0, iwrap=1,
 /
EOF

# ========================================================
# 7. Production
#    HMR, 4 fs timestep, 1 microsecond target
#    1 us = 1,000,000 ps = 250,000,000 steps at dt=0.004 ps
# ========================================================

cat <<EOF > 07_Prod_HMR_4fs.in
Production NPT with HMR and 4 fs timestep, 1 us target
 &cntrl
  imin=0, ntx=5, irest=1,
  ntc=2, ntf=2, tol=1.0e-7,
  dt=0.004, nstlim=250000000,
  ntt=3, gamma_ln=2.0, temp0=310.0, ig=-1,
  ntpr=125000, ntwr=250000, ntwx=125000,
  ntb=2, cut=10.0, ioutfm=1, ntxo=2,
  ntp=1, taup=1.0, barostat=2,
  ntr=0, iwrap=1,
 /
EOF

echo ">>> All 7 AMBER input files were created successfully."
