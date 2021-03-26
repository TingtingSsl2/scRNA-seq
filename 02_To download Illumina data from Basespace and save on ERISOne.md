# To download Illumina data from Basespace and save on ERISOne

- Sequencing data from Illumina is usually save in their Basespace platform. To download the data, a Basespace account is required (simple email registration).
- After log in Basespace, for each run, the run number can be found at the web address. This number is useful, copy and save.

To download the above data on ERISOne, follow the following code:
```
cd 
mkdir $HOME/bin
mkdir 2_data_Justin
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O $HOME/bin/bs
chmod u+x $HOME/bin/bs
wget https://api.bintray.com/content/basespace/BaseSpace-Copy-BIN/\$latest/linux/bscp?bt_package=develop -O $HOME/bin/bs-cp
bs auth # followed by pasting and adding address to web browser
bs download run -i 203846651 -o /PHShome/tz949/2_data_Justin/ # This is where that useful run number is used. 
```
