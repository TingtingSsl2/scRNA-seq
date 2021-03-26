# To install Cellranger on ERISOne

Prerequisites：
- Open an ERISOne account.
- If on BWH campus, internet should under Wifi: phswifi3.
- If WFH, set up VPN to connect to ERISOne. VPN set up see tutorial "BWH VPN set up".

Once log onto ERISOne, follow the following codes. 
```
cd
mkdir yard
cd yard
mkdir apps
cd apps
curl -o cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1616375849&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTYzNzU4NDl9fX1dfQ__&Signature=XlsdSdwGmuKJZZm1CudfNfVUF0hBJ7AnG63UK4I1rZQMdIem9ANlXvcuozJYGi1R9JMHEO8siNpqVHlUj9u9gyLx~DI4~XEYpheH6offBuY14jd7IFl4za1rXj6FG1eTLXHzDSe8bllLAgWmFp7DkZ4Nk0oq6OaFboDCOj150Iey7-~rOy07LfLljlgUD1SkFZRjBgZqkYz4L5WEeTZ-Kk7ZB6s~PSC75LyuX5gfHu8cyHj8uUBnXHuhZPrTmoyry4X7f5XRWjFdJeml62mRt2VhcX5bA7pA9rtkRxsywj7UBEUcGNt-VYhDiBbtbndFZCcU3jYpTxklgkwJTG2xJQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -zxvf cellranger-6.0.0.tar.gz
cd cellranger-6.0.0
cd /PHShome/tz949/yard/apps/cellranger-6.0.0
export PATH=/PHShome/tz949/yard/apps/cellranger-6.0.0:$PATH
which cellranger
```

After each log off and log in, cellranger needs to be added to the PATH as follows:
```
export PATH=/PHShome/tz949/yard/apps/cellranger-6.0.0:$PATH
```
