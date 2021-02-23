# elegansbrainmap


This is the accompanying code for 

Brittin, C. A. , Cook, S. J., Hall, D.H., Emmons, S. W., Cohen. N. A multiscale brain map derived from whole-brain volumetric reconstructions. Nature (in press). [biorxiv](https://doi.org/10.1101/2020.05.24.112870)

Source TrakEM2 data: [zenodo](https://zenodo.org/record/4383277#.X-wK5tZOk-I) or [wormwiring.org](http://wormwiring.org/) data.


## In preparation for advanced online publication, documentation will be under construction 21-24th February 2021

## Requirements:
Please check requirements.txt for details.

## Installation
```bash
pip install requrements.txt
```

## MySQL data
We found it convenient to store individual instances of membrane and synaptic contacts in a MySQL database, which can be downloaded from [here](https://zenodo.org/record/4383277#.X-wK5tZOk-I). MySQL offers several benefits, e.g., easier and faster to query data and less RAM intensive. Installing a MySQL database is straightforward (just do a simple web search for your particular OS). However, learning to navigate in SQL does have a bit of a learing curve. If you find that working with MySQL is a major hurdle, then reach out to us. If there is sufficient demand, then we can think about modifying the code to handle csv inputs.

### Import database
In your MySQL [create](https://www.digitalocean.com/community/tutorials/how-to-import-and-export-databases-in-mysql-or-mariadb) the databases 'N2U' and 'JSH'. [Download](https://zenodo.org/record/4383277#.X-wK5tZOk-I) databases.tar.bz2. Then import the databases as follows the adult:
```bash
cd databases
zcat adult_databases.sql.gz | mysql -u root -ppassword N2U
zcat l4_databases.sql.gz | mysql -u root -ppassword JSH
```
Substitute your user account with 'root' if appropriate and 'password' is the passord for the user account. 

### Accessing databases within the codebase
The codebase needs to 'login' everytime it accesses MySQL. Hence, your account info (i.e. username and password) needs to be physically accessible to the codebase. By default, we keep this info in '~/.my.cnf', but you can change the default file in the brainmap.db.connect module. For security reasons, we recommend keeping the .cnf file outside of your code repo (you do not want this info publicly available in your git repo). The .cnf file should have the following format
```
[client]
host=localhost
user=root-or-your-user-account-name
password=account-password
```
Once your MySQL databases and .cnf file are setup, you should be able run the scripts.

## Config file
The file configs/config.ini points to default mat and data file which are required by the scripts. Most scripts assume that the file configs/config.ini exists, however, for most scripts you can pass a different config file with the -c flag.  

## Generating paper figures
Most paper figures can be generated with the paper_figures.py wrapper. For example,
```bash
python paper_figures.py fig2a
python paper_figures.py ed1a
```
will generate Fig. 2a and Extended Data (ed) Fig. 1a, respectively. The file paper_figmap.csv maps figures to the appropriate scripts. 

## Preprocessing data
For convenience, we have provided preprocessed data in the data/ directory. In some instances where the preprocessed data is large (e.g. for Fig 2a), the figure generating the script will also do the necessary preproscessing. 

### Generate stochastic population model clusters
Default settings (e.g. Fig 1c)
```bash
python preprocess/cluster_population_models.py 
```
By default this will run 100 iterations in parallel across 10 cpus, to change the number of cpus and iterations
```bash
python preprocess/cluster_population_models.py -n 1 -i 1000
```
The default is σ=0.23, to change σ
```bash
python python preprocess/cluster_population_models.py -s 0.0,0.1,0.23,0.4
```
which will generate perturbations with σ=0.0, 0.1, 0.23 and 0.4. Use -h flag to see other adjustable settings. To plot the resulting cluster frequency matrices, take a look at analysis/cluster_population_plot_figure.py and scripts/fig1_cluster_results.py. 


