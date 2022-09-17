# Roadmap

* Half-yearly **releases**: **April and October**;
* Have a new documentation paper about every three to four years. The last one (ALF 2.0) was published on 2022.
* Development branches to be merged to master only if its **documentation** and CHANGELOG.md is up-to-date.


 To-do list for new releases:
 ---
  - check the pipeline tests pass
  - check the main documentation (paper) is up-to-date
  - check CHANGELOG is up-to-date
  - create ALF release branch from master
  - create corresponding pyALF branch
  - create corresponding Tutorial branch
  - update ALF website (links, and also 'news')


Next Goals
---

- Expand test pipelines to include simulations
- Improved overall usability:
  - Further develop pyALF
  - **Allow models to be input instead of coded**
- Time-dependent Hamiltonians
- Hybrid Monte Carlo with exact forces
- Interaction expansion (CT-INT)

**Check ideas in `BRAINSTORMING.md`**.


ALF 2.0
---
(_former "ALF 1.2"_)

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


