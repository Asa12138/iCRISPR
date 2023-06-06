#!/usr/bin/env python3
import os
add=os.popen("pwd")
a=add.read()
a_s=a.split('/')
smp=a_s[-2]
opt=open(smp+'.fa','a+')#‘a+’ ==a+r（可追加可读，文件若不存在就创建——潜在的所有 crispr
opt_cas=open(smp+'_cas.fa','a+')#从 opt 中抽出有cas注释的
existcrispr=os.popen("ls annotation_*")
e_s=existcrispr.read()#_s=string
e_l=e_s.split()#_l=list
e=[]#要遍历的含有 crispr 的 gff
for j in e_l:
    e.append(j[11:])

for k in e:
	ipt = os.popen("cat %s |cut -f 1,3,4,5,9|grep '_Crispr_\|spacer_'" %(k))#os.popen执行shell命令，结果返回到.txt#gff文件名要改#$k能传进去吗？os.popen("%s%"(k))
	i=ipt.readlines()
#	opt=open('norm_' + smp + k[:-4] + '.fa','w')#文件名要改

	cid=0

	for line in i:
	    if 'Crispr' in line:
	        top=line
	        cid+=1
	        sid=0
	        t=top.split()
	    else:
	        spacer=line
	        sid+=1
	        s=spacer.split()
	        last=s[-1]
	        l=last.split(';')
	        opt.write('>'+smp+'@'+str(s[0])+'@'+'CRISPR:'+str(cid)+'@'+str(t[2])+'-'+str(t[3])+'@'+str(l[1][5:])+'_'+str(sid)+'\n')
	        opt.write(l[0][9:]+'\n')

	ipt.close()
opt.close()#a+模式，指针在最后，read不到东西。所以要先close再open来read。


cas = os.popen("cat ../rawCas.fna|grep '>'|cut -d '|' -f 2|cut -d '_' -f 1|uniq") #输出含有 cas protein 的 array
cas_s=cas.read()
cas_l=cas_s.split('\n')[:-1]#最后一项是''（空字符串），不能包括在内

out=open(smp+'.fa','r')#a+模式，指针在最后，read不到东西。所以要先close再open来read。
crispr=out.readlines()

for x in cas_l:
	n=0
	while n<len(crispr):
		if x in crispr[n]:
			opt_cas.write(crispr[n]+crispr[n+1])
		n+=2

opt_cas.close()
out.close()
cas.close()#os.popen,open都要记得关！！！
