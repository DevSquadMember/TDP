# TDP4 - Factorisation LU et descente-remontée

## Equipe
- Dufrenne Alexis
- Saux Thomas

## Mise en route
Compiler les sources : `make build`

Pour lancer le programme en allouant dynamiquement des matrices avec des valeurs aléatoires :
`make NP=<nb_procs> run <size>`

Attention, la taille des matrices doit être proportionnelle au nombre de processeurs utilisés.

## Tests

Pour vérifier que l'algorithme trouve bien le vecteur solution X au problème :
`make NP=<nb_procs> check <size>`