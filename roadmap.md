
Main points from March meetings:

* Adopt half-yearly **releases**: **1st April and 1st October**;
* Have a new documentation paper about every three years. The next one should be online by June 2020.
* ALF's April 2020 release = ALF 2.0 (instead of placeholder "ALF 1.2").

Also:

* A development branch is to be merged to master only if its **documentation** and CHANGELOG.md is up-to-date.


 To-do list for new releases:
 ---
  - check the main documentation (paper) is up-to-date
  - check changelog is up-to-date
  - create ALF release branch from master
  - create corresponding pyALF branch
  - create corresponding Tutorial branch
  - update ALF website (links, and also 'news')


ALF 2.0
---

- Parallel tempering
- Global updates / global tau
- Projective approaches
- Continuous fields
- Langevin
- Maxent
- More models
- Trotter symmetric
- Predefined structures.
The last two still need some work from Fakher, and Langevin has to be merged in, but then the code will be essentially ready for release and we only have to worry about the documentation (paper).

ALF 3.0
---

We intend to have:
- Improved overall usability
- Improved I.O., HDF5
- Improved post-processing
- Implementation with classes
- Trotter options
- Time-dependent Hamiltonians
- Hybrid Monte Carlo with exact forces
- Rényi entropies
- Interaction expansion (CT-INT).
