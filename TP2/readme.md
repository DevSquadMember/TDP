# TDP2 - Gravitation universelle

## Equipe
- Dufrenne Alexis
- Saux Thomas

## Mise en route
- `mpirun -np 4 cmake-build-debug/main`
- `plot "res.dat" using 2:3:4:5 with vectors, "res.dat" using 6:7:8:9 with vectors`

## Tâches
- communications dans l'anneau avec double-buffer (un en écriture et envoi, l'autre en réception)
- structure des données (structure de la masse, et données à envoyer : m et p)
- calcul des forces, vitesses, nouvelle position
- affichage, rendu + fichier .dot

## A FAIRE
- make/cmake
- commentaires
- sortie graphique en gnuplot/pythons
- entrée : nombre d'itérations en temps
- fichier :
n (fois)
m px py vx vy
........
- code de test qui compare le séquentiel en parallèle

## Performances
- fixer n et augmenter P
- utiliser n assez grand, n > 5000
- être clair dans les conditions expérimentales (nb noeuds, nb process/noeud)
- imposer des conditions d'entrée : n=alpha*p par exemple

Pour chaque itération t -> t + dt
  calcul des forces (p tours)
  calcul du dt global
  mise à jour des particules avec dt global
