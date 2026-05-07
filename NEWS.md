# DeadWood_DWDDecay News

## v0.0.1

- Initial release.
- Receives fallen snags from `DeadWood_snagDecay`, remaps decay class to the
  DWD scale, and advances DWD annually via a Markov transition matrix.
- Stochastically removes fully decomposed DWD records each timestep.
