# Roadmap


Merges
---

* Development branches to be merged to master only when their 
  - **documentation** reflects all changes
  - **CHANGELOG.md** mentions all breaking changes and critical bug fixes


Release cycle
---

* Half-yearly releases: **April and October**
* Main release (new documentation paper): ~**every three to four years**

Last main release: ALF 2.0, published on 2022 (2020 on ArXiV).


 To-do list for new releases:
 ---
  - check the pipeline tests pass
  - check the main documentation (paper) is content-wise up-to-date (all changes are reflected in the text)
  - check CHANGELOG is up-to-date
  - create ALF release branch from master, then
    - add version number to header of README (ALF -> ALF #.#)
    - update main documentation (paper — LaTeX variables such `\ALFver` correctly set)
  - create corresponding pyALF branch (or add remark to its README it's compatible with the new version)
  - create corresponding Tutorial branch (or add remark to its README it's compatible with the new version)
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


