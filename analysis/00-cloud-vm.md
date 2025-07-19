# Run scripts on DO

Setup DO server (ubuntu 24.04 LTE, 8 vCPU 16G RAM):

```bash
export IP=142.93.151.236
ssh root@$IP

curl -sSL https://repos.insights.digitalocean.com/install.sh | sudo bash

# https://github.com/eddelbuettel/r-ci/blob/master/docs/run.sh
curl -OLs https://eddelbuettel.github.io/r-ci/run.sh && chmod 0755 run.sh
./run.sh bootstrap

R -q -e "install.packages(c('bSims', 'detect', 'qs', 'rconfig', 'extraDistr'))"

exit
```

On local machine, push files up to server:

```bash
rsync -azP \
    --exclude '.git' \
    /Users/Peter/git/github.com/psolymos/bsims-tests \
    root@$IP:/root

rsync -azP \
    /Users/Peter/git/github.com/psolymos/bsims-tests/sqpad-paper/analysis \
    root@$IP:/root

```

Run R script:

```bash
ssh root@$IP

# tmux new -s mysession
# tmux kill-ses -t mysession
# tmux ls
# tmux a -t mysession
# Ctrl+b then d - detach
# Ctrl+b then x - kill


cd bsims-tests/sqpad-paper/analysis

Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 1 # running on 0
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 2 # running on 1
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 3 # running on 2
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 4 # running on 3

Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 5 # running on 4
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 6 # running on 5
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 7 # running on 6
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 8 # running on 7

Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 9 # running on 8
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 10 # running on 9
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 11 # running on 10
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 12 # running on 11

Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 13 # running on 12
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 14 # running on 13
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 15 # running on 14
Rscript --vanilla 03-bsims-estimation.R --B 2000 --notify --set 16 # running on 15


# this is the current
ls -al /root/bsims-tests/sqpad-paper/analysis/_tmp/est_conv_mc
```

Push notifications:

- install app, see: <https://docs.ntfy.sh/>
- set up a topic, e.g. `a8m_cht_alerts`
- send alert `curl -d "Run successful" ntfy.sh/a8m_cht_alerts`

In R:

```R
topic <- "a8m_cht_alerts"
msg <- paste("TEST Finished", YEAR, "@", .POSIXct(Sys.time(), "America/Edmonton"))
system2("curl", c("-d", sprintf("\"%s\"", msg), sprintf("ntfy.sh/%s", topic)))
```

Once finished, copy results back (keep timestamps too):

```bash
# this is the current
rsync -rt \
    root@$IP:/root/bsims-tests/sqpad-paper/analysis/_tmp/est_conv_mc \
    /Users/Peter/git/github.com/psolymos/bsims-tests/sqpad-paper/analysis/_tmp

rsync -rt \
    root@$IP:/root/analysis/_tmp/paired_mc \
    /Users/Peter/git/github.com/psolymos/bsims-tests/sqpad-paper/analysis/_tmp

```
