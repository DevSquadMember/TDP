# TDP3 - Produit matriciel Fox

## Equipe
- Dufrenne Alexis
- Saux Thomas

## Mise en route
Compiler les sources : `make build`

Le programme charge les fichiers `a.txt` et `b.txt` pour les matrices A et B respectivement.

Pour lancer le programme depuis les fichiers :
`make NP=<nb_procs> run`

Il est possible de générer de manière aléatoires des matrices avec des valeurs comprises entre -100 000 et 100 000
en précision 10e-6, pour ce faire :
`make generate <size>`

Il est également possible d'exécuter le programme sans passer par des fichiers mais en allouant dynamiquement des matrices
avec des valeurs aléatoires :
`make NP=<nb_procs> run <size>`