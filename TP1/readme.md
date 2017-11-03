# TDP1 - Mise en oeuvre des BLAS

## Equipe
- Dufrenne Alexis
- Saux Thomas

## Mise en route
- `make` pour compiler les sources
- `./driver` pour lancer le programme

## Choix des méthodes à lancer
Dans le fichier source `driver.c`, plusieurs méthodes ont été implémentées afin de pouvoir tester les performances, comme par exemple pour ddot ou dgemm.
Il suffit dans le corps du `main` de commenter ou décommenter les appels aux fonctions.
L'exécution du programme `driver` affiche les résultats et sauvegarde dans `res.dot` des valeurs exploitables par `gnuplot`