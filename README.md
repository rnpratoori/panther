PANTHER
=====

**P**hase-separation **AN**alysis **T**ool for **H**eterogeneous polym**ER**s (PANTHER) is a MOOSE-based application for studying phase separation behavior in block copolymers and crosslinked polymers.

# Compilation instructions:
1. Install MOOSE following these instructions: https://mooseframework.inl.gov/getting_started/installation/index.html
2. Clone this repository:
   ```
   git clone https://github.com/rnpratoori/panther
   ```
4. Compile the repository:
   ```
   cd panther
   make -j 6
   ```
5. Test if the repository is working fine:
   ```
   ./run_tests
   ```
