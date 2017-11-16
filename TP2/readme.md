# TDP2 - Gravitation universelle

## Equipe
- Dufrenne Alexis
- Saux Thomas

## Mise en route
Compiler les sources :
`make`

Lancer la version séquentielle :
`make run_seq <nb_iterations> <fichier planètes> <générer graphe : 0 ou 1>`

Lancer la version parallèle :
`make NP="<nb_proc>" run_par <nb_iterations> <fichier planètes> <générer graphe : 0 ou 1>`

Le nombre de particules/planètes doit être un multiple du nombre de processeurs.

Générer un fichier de particules aléatoires :
`make generate <nb_particules> <nom du fichier de sortie>`

# Test des performances
Un programme a été développé afin de pouvoir mesurer les différences de performances 
entre la simulation en version séquentielle et celle en version parallèle.

La configuration des tests se fait dans le fichier `compare.c`, où on peut modifier :
- le nombre de processeurs max utilisés : `NB_PROC_MAX`
- le nombre de particules : `NB_PARTICULES`
- le nombre d'itérations de départ : `NB_MIN_ITERATIONS`
- le nombre max d'itérations : `NB_MAX_ITERATIONS`
- le coefficient à appliquer pour incrémenter le nombre d'itérations : `COEF`
