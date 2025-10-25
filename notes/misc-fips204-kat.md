
## KAT for FIPS-203, FIPS-204, and SP800-208

This repository contains a set of randomly generated Known Answer Test (KAT) vectors 
for FIPS-203 and FIPS-204: https://github.com/post-quantum-cryptography/KAT


## Steps to Prepare KAT File
```bash
$ pwd
/Users/dprasad/projects/dicp-local/KAT/MLDSA

$ cp kat_MLDSA_44_hedged_raw.rsp ../../ml-dsa/kats/misc-keygen44-kats.txt
$ cp kat_MLDSA_65_hedged_raw.rsp ../../ml-dsa/kats/misc-keygen65-kats.txt
$ cp kat_MLDSA_87_hedged_raw.rsp ../../ml-dsa/kats/misc-keygen87-kats.txt

```