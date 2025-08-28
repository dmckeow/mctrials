# mctrials
Practice for manual curation of TEs in Drosophila
## Workflow

### Setup
I was working temporarily in /home/csic/gcy/jgp/extra_storage/dean/mctrials/mctrials

```{bash}
conda create -n bioinf -y
conda activate bioinf
conda install -c bioconda seqkit
```


### Running TEAid through MCHelper
* Need to run TEAid independently through MCHelper - but only the version forked by Adrian
* The library input MUST be the NR curated library output from MCHelper


### Getting data into TEammo
I ran MCHelper on pregenerated raw RM2 libraries, using the CESGA cluster - now I will transfer the results to a local workstation

```{bash}
rsync -avz csgcyjgp@ft3.cesga.es:/mnt/netapp2/Store_csgcyjgp/dean/mctrials/mctrials/data/ mctrials_data/
```

Then we linked the data to the data folder where TEammo looks and stores outputs
* The toy data is linked with `ln -s ../teammo_data/ data`
```{bash}
cd ~/TEammo
rm -f data
ln -s ../mctrials_data/ data
```


### Running TEammo
* On a local workstation, ssh through own computer
* Login and installation setup by Adrian, but likely it was the conda setup for TEammo he did
  
```{bash}
ssh -L 3001:localhost:3001 dmckeow@161.111.135.77

cd TEammo
# If you just run the Shiny app in terminal without redirect it is unstable
Rscript TEammo_app.R > teammo.log 2>&1 &

# In case a process was not cancelled properly and it is stilling running
lsof -i :3001
kill -9 <pid>
```