### Syncing
To set up use `rclone config`, then it's pretty easy and login through broswer.

Ideally use crontab, but that ain't working on CX3.

For RDS.
```
rclone sync /rds/general/user/jl24018/home/projects/openmm-md/data/241010_FoldingUponBinding/ onedrive:/data/241010_FoldingUponBinding/
```