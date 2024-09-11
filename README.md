**💻🧬 BLAST 2**\
_Basic Local Alignment Search Tool 2_
==============
<ins>Étudiante</ins> : Louise LAM   

## Pour commencer
Ces instructions vont vous permettre d'obtenir une version fonctionnelle de mon projet court Blast2 sur votre machine locale. Lisez les instructions suivantes pour une expérience optimale.

## Prérequis
Ce projet a été testé avec `Python 3.12.5` et un environnement virtuel `conda 24.7.1`. L’installation de ces versions sont recommandées pour pouvoir lancer le projet correctement.

## Installation
1. Cloner le dépot du projet depuis le répertoire parent. Dans une console UNIX, utiliser la commande `git clone`.

2. Se placer dans le répertoire `blast2` : 
```
cd blast2
```

3. Créer l'environnement virtuel `conda` à partir du fichier `environment.yml` à l'aide de la commande :
```
conda env create -f environment.yml
```



## Utilisation
1. Activer l'environnement virtuel `conda` à l'aide de la commande : 
```
conda activate blast2-env
```
2. Lancer blast2 sur la séquence de la protéine P53 humaine :
```
python src/main.py -d data/subset.fasta -q data/p53_human.fasta
```

