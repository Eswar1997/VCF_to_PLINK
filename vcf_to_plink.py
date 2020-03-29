#Import required packages
import sys,getopt
import os
import pandas as pd
import numpy as np
import re
import time

##################################Map file conversion fumction########################################
def vcf_to_map(vcf_file,output_file_name):
    map_file=pd.DataFrame(columns=["CHR","BP","ID"])
    ext_file=pd.DataFrame(columns=["CHR","BP","ID","REF","ALT"])
    map_file=vcf_file.iloc[:,0:3].copy()
    #Uncomment the following line if you want the allele information in a .ext file
    #ext_file=vcf_file.iloc[:,0:5].copy()
    map_file['DIST'] = int('0')
    map_file.columns = ['CHR','BP','ID','DIST']
    map_file = map_file.reindex(columns=['CHR','ID','DIST','BP'])
    map_file = map_file.replace(to_replace ='chr', value = '', regex=True)
    map_file.to_csv(output_file_name+'.map',sep="\t",index=False,header=False)
    #Uncomment the following line for .ext file
    #ext_file.to_csv(output_file_name+'.ext',sep="\t",index=False,header=True)
    print("First fewlines of mapfile is:")
    print(map_file.head())
    print("Map file created successfully")

############################################Ped file Conversion###########################################
def vcf_to_ped(vcf_file,output_file_name):
    vcf_file_reqd = vcf_file.iloc[:,9:].copy()
    vcf_file_trans = vcf_file_reqd.transpose()
    columns=len(vcf_file_trans.columns)
    rows=vcf_file_trans.shape[0]
    #print(columns)
    #print(rows)
    ped_file=vcf_file_trans.copy()
    cols=ped_file.columns
    prev_len=len(cols)
    print("Number of columns in pedfile should be 6+2*"+str(prev_len)+"=",(6+(2*prev_len)))
    ped_file.index.name = 'FID'
    ped_file.reset_index(inplace=True)
    ped_file['IID']=ped_file['FID']
    ped_file['PID']=ped_file['MID']=ped_file['Sex']=str(0)
    ped_file['Pheno']=str(-9)
    cols = list(ped_file)
    cols.insert(1, cols.pop(cols.index('IID')))
    cols.insert(2, cols.pop(cols.index('PID')))
    cols.insert(3, cols.pop(cols.index('MID')))
    cols.insert(4, cols.pop(cols.index('Sex')))
    cols.insert(5, cols.pop(cols.index('Pheno')))
    ped_file = ped_file.loc[:, cols]
    #print(ped_file.head())
    ###########################Convering genotypes according to plink formats##############################
    ped_file=ped_file.replace({r'^\.\/\.*':'0/1',r'^1/1*':'2/2',r'^0/0*':'0/0',r'^0/1*':'1/0'},regex = True)
    ped_file=ped_file.astype(dtype=str, copy=True, errors='raise')
    cols=ped_file.columns
    ped_file= pd.concat([ped_file[x].str.split('/', expand=True) for x in cols], axis=1, keys=ped_file.columns)
    ped_file.columns = ped_file.columns.map(lambda x: '_'.join((str(x[0]), str(x[1]))))
    #ped_file = ped_file.replace({r':*':''},regex=True)
    #pedfile_updated=ped_file.replace(to_replace =r'\1:*', value = '', regex = True)
    #ped_file_updated = ped_file(lambda x:':'.join(str(x[0]).split(':')[:-1]))
    ped_file=ped_file.replace({r'^0:*': 0,r'^1:*':1,r'^2:*':2}, regex = True)
    #ped_file=ped_file.replace(to_replace =r'^1:*', value = 1, regex = True)
    #ped_file=ped_file.replace(to_replace =r'^2:*', value = 2, regex = True)
    #print(ped_file.head())
    print("Total columns in ped files:" ,ped_file.shape[1])
    ped_file.to_csv(output_file_name+'.ped',sep="\t",index=False,header=False)

#######################The convert vcf to ped and map controlling function#############################
def convert_vcf(input_file_name,output_file_name):
    #print("VCF file is :" +input_file_name)
    #print("Plink formatted file is :" +output_file_name)
    pattern = re.compile("^#CHROM")
    for i, line in enumerate(open(input_file_name)):
        for match in re.finditer(pattern, line):
            #Uncomment the following line if you want to know where the vcf header is ending
            print ('VCF headers till line %s: %s' %(i+1, match.group()))
            start_at=int(i)

    vcf_file = pd.read_csv(input_file_name,header = start_at,sep='\s',na_filter = True,skip_blank_lines= True, engine='python')
    #print(vcf_file.head())
    print("Number of markers in the given vcf file is ",vcf_file.shape[0])
    print("Number of samples in the given vcf file is ",vcf_file.shape[1])

    pre= time.perf_counter()
    vcf_to_map(vcf_file,output_file_name)
    post= time.perf_counter()
    print("Time to generate mapfile:"+str(post-pre)+"s")
    pre= time.perf_counter()
    vcf_to_ped(vcf_file,output_file_name)
    post= time.perf_counter()
    print("Time to generate pedfile:"+str(post-pre)+"s")


 ##############################Python main function - getting command line input#############################

def main(argv):
   print("VCF to PLINK format converter")
   print("Written by Alagu Sankareswaran. Suggestions are welcome.")
   print("EmailID: sankareswaran@ccmb.res.in / sankareswaransubramanian@gmail.com")
   print("Provide only file name. Extension will be taken care automatically")
   input_file_name = ''
   output_file_name= ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print ('vcf_to_plink.py -i <vcf_file_name> -o <plink_format_file_name>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('vcf_to_plink.py -i <vcf_file_name> -o <plink_format_file_name>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         input_file_name = arg
      elif opt in ("-o","--ofile"):
          output_file_name = arg

   print('VCF file name is "', input_file_name)
   print('Plink file name is "', output_file_name)
   input_file_name = input_file_name +".vcf"
   convert_vcf(input_file_name,output_file_name)

if __name__ == "__main__":
   main(sys.argv[1:])
