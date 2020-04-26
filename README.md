# VCF_TO_PLINK

The __vcf_to_plink.py__ script will convert Unphased vcf file to plink format and it is written in _Python3_. The converted file can be used in plink with --file option.


## Usage of the script

The program can be run using the following command

```python vcf_to_plink.py -i <vcf_file_name> -o <plink_format_file_name>```

To get the above help provide -h in the terminal

```python vcf_to_plink.py -h```

Two file names are required for the script to run successfully. Input file succeding '-i' is the vcf file name you want to convert. '-o' is the plink formatted output file name. Please don't add any extension after the file name such as '.vcf' or '.ped/.map'. These are all taken care by the script. 

