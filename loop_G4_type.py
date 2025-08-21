import sys
import pandas as pd

loop_file=sys.argv[1]
G4_anchor_file=sys.argv[2]
TDP43_anchor_file=sys.argv[3]
anchor_promoter_file=sys.argv[4]
promoter_anchor_file=sys.argv[5]
enhancer_anchor_file=sys.argv[6]
loop_type_file=sys.argv[7]

G4_anchor=[]
f=open(G4_anchor_file)
for line in f:
	line=line.strip("\n")
	array=line.split("\t")
	G4_anchor.append(array)
f.close()

TDP43_anchor=[]
f=open(TDP43_anchor_file)
for line in f:
	line=line.strip("\n")
	array=line.split("\t")
	TDP43_anchor.append(array)
f.close()


promoter_anchor=[]
f=open(promoter_anchor_file)
for line in f:
	line=line.strip("\n")
	array=line.split("\t")
	promoter_anchor.append(array)
f.close()

enhancer_anchor=[]
f=open(enhancer_anchor_file)
for line in f:
	line=line.strip("\n")
	array=line.split("\t")
	enhancer_anchor.append(array)
f.close()


anchor_promoter={}
df=pd.read_csv(anchor_promoter_file,sep="\t",header=None)
print(df)


f=open(loop_file)
out=open(loop_type_file,'w')
for line in f:
	line=line.strip("\n")
	array=line.split("\t")
	temp=[array[0],array[1],array[2]]
	type_list=[]
	if temp in G4_anchor:
		type1="G4"
	else:
		type1="None"
	type_list.append(type1)
	type3=""
	if temp in TDP43_anchor:
		type3="MAX"
	else:
		type3="None"
	type5=""
	if temp in promoter_anchor:
		type5="P"
	else:
		type5="None"
	type7=""
	if temp in enhancer_anchor:
		type7="E"
	else:
		type7="None"
	type2=""
	temp=[array[3],array[4],array[5]]
	if temp in G4_anchor:
		type2="G4"
	else:
		type2="None"
	type4=""
	if temp in TDP43_anchor:
		type4="MAX"
	else:
		type4="None"
	type6=""
	if temp in promoter_anchor:
		type6="P"
	else:
		type6="None"
	type8=""
	if temp in enhancer_anchor:
		type8="E"
	else:
		type8="None"
	df_anchor1=df[(df[10]==array[0])&(df[11]==int(array[1]))&(df[12]==int(array[2]))]
	df_anchor2=df[(df[10]==array[3])&(df[11]==int(array[4]))&(df[12]==int(array[5]))]
	promoters1=""
	promoters2=""
	if len(df_anchor1)==0:
		promoters1=""
	else:
		promoters1_list=df_anchor1[4].unique()
		for i in promoters1_list:
			promoters1=promoters1+i+","
	if len(df_anchor2)==0:
		promoters2=""
	else:
		promoters2_list=df_anchor2[4].unique()
		for i in promoters2_list:
			promoters2=promoters2+i+","
	out.write(line+"\t"+type1+"\t"+type2+"\t"+type3+"\t"+type4+"\t"+type5+"\t"+type6+"\t"+type7+"\t"+type8+"\t"+promoters1+"\t"+promoters2+"\n")
out.close()
f.close()
