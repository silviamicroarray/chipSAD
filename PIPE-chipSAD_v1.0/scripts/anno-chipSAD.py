import os
import sys
from optparse import OptionParser
from math import sqrt
from operator import itemgetter
from scipy import stats
import fnmatch 

def main():
    opts=parse_cmd_line()
    
    table(opts) 
    
    
    
def parse_cmd_line():
    #Parse the command line and return a parameters structure.

    usage = "usage: %prog [options] [structure]"
       
    parser = OptionParser(usage)
    parser.add_option("-o", "--outdir", action="store", dest="outdir", default=".", type="string", help="directory for output file. Default is current directory.")
    parser.add_option("-u", "--thup", action="store", dest="thup", default=1.00, type="float", help="up threshold for DE transcripts. Default is 1.00.")
    parser.add_option("-d", "--thdown", action="store", dest="thdown", default=-1.00, type="float", help="down threshold for DE tramscripts. Default is -1.00.")
    


    (options, args) = parser.parse_args()
    if len(args) != 4:
        parser.error("The chipSAD results file for each strand, a probe file, a .gbk file are required.")

    if not os.path.isdir(options.outdir):
        parser.error("Not a valid directory: '%s'." % options.outdir)

    options.entity1 = args[0]
    options.entity2= args[1]
    options.entity3=args[2]
    options.entity4=args[3]
   
    
    return options
    
    
    
def table(options):
    file_in_name1 = os.path.join(options.outdir, options.entity1)
    file_in_name2 = os.path.join(options.outdir, options.entity2)
    file_probes_name1 = os.path.join(options.outdir, options.entity3)
    file_gbk_name = os.path.join(options.outdir, options.entity4)
    

    file_out_base = os.path.splitext(os.path.basename(options.entity1))[0]
    file_out_base2 = os.path.splitext(os.path.basename(options.entity2))[0]
    
    file_out_name=os.path.join(options.outdir,'tot_results_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_out_name1f = os.path.join(options.outdir,'F_UTR_results_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_out_name2f = os.path.join(options.outdir,'F_UTR_over_results_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_out_name3f = os.path.join(options.outdir,'F_antisense_results_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_out_name4f = os.path.join(options.outdir,'F_operons_results_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_out_name5f = os.path.join(options.outdir,'F_intergenic_results_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_artemis_utr_namef = os.path.join(options.outdir,'F_UTR_artemis_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_artemis_utr_over_namef = os.path.join(options.outdir,'F_UTR_over_artemis_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_artemis_antisense_namef = os.path.join(options.outdir,'F_antisense_artemis_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_artemis_operons_namef = os.path.join(options.outdir,'F_operons_artemis_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    file_artemis_inter_namef = os.path.join(options.outdir,'F_intergenic_artemis_%s_th_u%s_d%s.txt' % (file_out_base, options.thup, options.thdown))
    
    file_out_name1r = os.path.join(options.outdir,'R_UTR_results_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_out_name2r = os.path.join(options.outdir,'R_UTR_over_results_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_out_name3r = os.path.join(options.outdir,'R_antisense_results_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_out_name4r = os.path.join(options.outdir,'R_operons_results_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_out_name5r = os.path.join(options.outdir,'R_intergenic_results_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_artemis_utr_namer = os.path.join(options.outdir,'R_UTR_artemis_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_artemis_utr_over_namer = os.path.join(options.outdir,'R_UTR_over_artemis_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_artemis_antisense_namer = os.path.join(options.outdir,'R_antisense_artemis_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_artemis_operons_namer = os.path.join(options.outdir,'R_operons_artemis_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    file_artemis_inter_namer = os.path.join(options.outdir,'R_intergenic_artemis_%s_th_u%s_d%s.txt' % (file_out_base2, options.thup, options.thdown))
    
    file_in1 = file(file_in_name1, 'r')
    file_in2 = file(file_in_name2, 'r')
    file_probes = file(file_probes_name1, 'r')
    file_gbk= file(file_gbk_name, 'r')
    
    file_out = file(file_out_name, 'w')
    file_utrf = file(file_out_name1f, 'w')
    file_utr_overf = file(file_out_name2f, 'w')
    file_antisensef = file(file_out_name3f, 'w')
    file_operonsf = file(file_out_name4f, 'w')
    file_interf = file(file_out_name5f, 'w')
    file_artemis_utrf = file(file_artemis_utr_namef, 'w')
    file_artemis_utr_overf = file(file_artemis_utr_over_namef, 'w')
    file_artemis_antisensef = file(file_artemis_antisense_namef, 'w')
    file_artemis_operonsf = file(file_artemis_operons_namef, 'w')
    file_artemis_interf = file(file_artemis_inter_namef, 'w')
    
    file_utrr = file(file_out_name1r, 'w')
    file_utr_overr = file(file_out_name2r, 'w')
    file_antisenser = file(file_out_name3r, 'w')
    file_operonsr = file(file_out_name4r, 'w')
    file_interr = file(file_out_name5r, 'w')
    file_artemis_utrr = file(file_artemis_utr_namer, 'w')
    file_artemis_utr_overr = file(file_artemis_utr_over_namer, 'w')
    file_artemis_antisenser = file(file_artemis_antisense_namer, 'w')
    file_artemis_operonsr = file(file_artemis_operons_namer, 'w')
    file_artemis_interr = file(file_artemis_inter_namer, 'w')
    
    
    thup=options.thup
    thdown=options.thdown
    
    
    print "gbk file elaboration"
    
    #parsing gbk
    gbk_array=[]
    start_array_gbk_r=[]
    end_array_gbk_r=[]
    start_array_gbk_f=[]
    end_array_gbk_f=[]
    name_array_gbk_r=[]
    name_array_gbk_f=[]
    start_array_gbk_glob=[]
    end_array_gbk_glob=[]
    name_array_gbk_glob=[]
    #location_array_gbk=[]
    #synonymus_array_gbk=[]
    strands=''
    strand_glob_array=[]
    for line in file_gbk:
        
        if line.startswith('     gene'):
            location=line.split('            ')[1].split('\n ')[0]            
            
            if fnmatch.fnmatch(location, '*complement*'):
                #loc=location.partition('complement')[2].split('(')[1].split(')')[0].split('\n')[0]
                start=int(location.split('(')[1].split('..')[0])
                end=int(location.split('(')[1].split('..')[1].split(')')[0])
                #location_array_gbk.append(loc)
                start_array_gbk_r.append(start)
                end_array_gbk_r.append(end)
                start_array_gbk_glob.append(start)
                end_array_gbk_glob.append(end)
                #file_out.write("gbk: start %s, end %s\n" %(start, end))
                #file_out.write("gbk: c  location %s\n" %(loc))
                #file_out.write("%s\n" %(loc))
                strands='R'
            else:
                #loc=location.split('\n')[0]
                #location_array_gbk.append(loc)
                start=int(location.split('..')[0])
                end=int(location.split('..')[1].split('\n')[0])
                start_array_gbk_f.append(start)
                end_array_gbk_f.append(end)
                start_array_gbk_glob.append(start)
                end_array_gbk_glob.append(end)
                #file_out.write("gbk: start %s, end %s\n" %(start, end))
                #file_out.write("gbk:   location %s\n" %(loc))
                #file_out.write("%s\n" %(loc))
                strands='F'
                
            #print start
            #print end
        if line.startswith('                     /locus_tag'):
            
            name=(line.split('="')[1].split('"')[0])
            if strands=='F':
                name_array_gbk_f.append(name)
                name_array_gbk_glob.append(name)
                strand_glob_array.append('>')
            if strands=='R':
                name_array_gbk_r.append(name)
                name_array_gbk_glob.append(name)
                strand_glob_array.append('<')
            strands=''
            
                
    
    
    #print ("gene F %s\n" %(len(start_array_gbk_f)))
    #print ("gene R %s\n" %(len(start_array_gbk_r)))
    #print ("name F %s\n" %(len(name_array_gbk_f)))
    #print ("name R %s\n" %(len(name_array_gbk_r)))
    
    print "data elaboration"
    #elaborazione dati
    max_array1=[]
    #average_array=[]
    pvalue_array1=[]
    #file_array=[]
    start_array1=[]
    end_array1=[]
    line_array1=[]
    ttest_array1=[]
    for line in file_in1:
            if (line.startswith('***')):
                sent=line.split()
                #print sent
                kconts= [sent.count(w) for w in sent]
                #file_out.write("count %d\n" %(kcont))
                kcont=len(kconts)
                #print kcont
                if kcont==14:
                    st=int(line.split()[2])
                    end=int(line.split()[4])
                    start_array1.append(st)
                    end_array1.append(end)
                    line_array1.append(line)
                    maxv=float(line.split()[11]) #non sul max ma sul valore medio
                    ttest=float(line.split()[12])
                    pval=float(line.split()[13])
                    max_array1.append(maxv)
                    ttest_array1.append(ttest)
                    pvalue_array1.append(pval)
                    
    max_array2=[]       
    start_array2=[]
    end_array2=[]
    line_array2=[]
    pvalue_array2=[]
    ttest_array2=[]
    for line in file_in2:
            if (line.startswith('***')):
                sent=line.split()
                #print sent
                kconts= [sent.count(w) for w in sent]
                #file_out.write("count %d\n" %(kcont))
                kcont=len(kconts)
                #print kcont
                if kcont==14:
                    st=int(line.split()[2])
                    end=int(line.split()[4])
                    maxv=float(line.split()[11]) #non sul max ma sul valore medio
                    ttest=float(line.split()[12])
                    pval=float(line.split()[13])
                    max_array2.append(maxv)
                    start_array2.append(st)
                    end_array2.append(end)
                    line_array2.append(line)
                    pvalue_array2.append(pval)
                    ttest_array2.append(ttest)
            
    
    
    
    start_probe_arrayf=[]
    end_probe_arrayf=[]
    status_arrayf=[]
    nome_probe_arrayf=[]
    value_arrayf=[]
    start_probe_arrayr=[]
    end_probe_arrayr=[]
    status_arrayr=[]
    nome_probe_arrayr=[]
    value_arrayr=[]
    kcont=0
    
    for line in file_probes:
        
            if (line.startswith('F')):
                sent=line.split()
                kconts= [sent.count(w) for w in sent]
                #file_out.write("count %d\n" %(kcont))
                kcont=len(kconts)
                if kcont==8:
                    stp=int(line.split()[3])
                    endp=int(line.split()[4])
                    status=line.split()[1]
                    nome=line.split()[5]
                    val=float(line.split()[6])
                    
                    start_probe_arrayf.append(stp)
                    end_probe_arrayf.append(endp)
                    status_arrayf.append(status)
                    nome_probe_arrayf.append(nome)
                    value_arrayf.append(val)
                   
                if kcont==7:
                    stp=int(line.split()[2])
                    endp=int(line.split()[3])
                    status=line.split()[1]
                    nome=line.split()[4]
                    val=float(line.split()[5])
            
                    start_probe_arrayf.append(stp)
                    end_probe_arrayf.append(endp)
                    status_arrayf.append(status)
                    nome_probe_arrayf.append(nome)
                    value_arrayf.append(val)
                    
    
        
            if (line.startswith('R')):
                sent=line.split()
                kconts= [sent.count(w) for w in sent]
                #file_out.write("count %d\n" %(kcont))
                kcont=len(kconts)
                if kcont==8:
                    stp=int(line.split()[3])
                    endp=int(line.split()[4])
                    status=line.split()[1]
                    nome=line.split()[5]
                    val=float(line.split()[6])
                    
                    start_probe_arrayr.append(stp)
                    end_probe_arrayr.append(endp)
                    status_arrayr.append(status)
                    nome_probe_arrayr.append(nome)
                    value_arrayr.append(val) 
                   

                if kcont==7:
                    stp=int(line.split()[2])
                    endp=int(line.split()[3])
                    status=line.split()[1]
                    nome=line.split()[4]
                    val=float(line.split()[5])
                    
                    start_probe_arrayr.append(stp)
                    end_probe_arrayr.append(endp)
                    status_arrayr.append(status)
                    nome_probe_arrayr.append(nome)
                    value_arrayr.append(val)  
                    
    #print len(start_array1)
    #print len(start_array2)
    #print len(start_probe_arrayf)
    #print len(start_probe_arrayr)
    
    #print start_array_gbk_f
    
    #print start_probe_arrayf
    #for z in range(len(status_arrayf)):
       # file_utr.write("%s" %(status_arrayf[z]))
    start_utr=[]
    end_utr=[]
    max_utr=[]
    start_utr_over=[]
    end_utr_over=[]
    max_utr_over=[]
    start_anti=[]
    end_anti=[]
    max_anti=[]
    start_operonsf=[]
    end_operonsf=[]
    max_operonsf=[]
    start_interf=[]
    end_interf=[]
    max_interf=[]
    temp_array_ref=[]
    temp_array_refr=[]
    temp_array_ref1=[]
    temp_array_ref0=[]
    temp_glob_st=[]
    temp_glob_end=[]
    temp_glob_strand_st=[]
    temp_glob_strand_end=[]
    temp_array_st=[]
    temp_array_end=[]
    orfs=''
    orfover=''
    countCf=0
    countCr=0
    countNf=0
    countNr=0
    countTf=0
    countTr=0
    countGenef=0
    countGener=0
    count_gene_f=0
    count_gene_r=0
    
    i=0
    j=0
    k=0
    kutr=0
    kutrover=0
    kanti=0
    kop=0
    kinter=0
    count_nobuchi=0
    for i in range(len(start_array1)):
        countCf=0
        countCr=0
        countNf=0
        countNr=0
        countTf=0
        countTr=0
        countGenef=0
        countGener=0
        count_gene_f=0
        count_gene_r=0
        count_nobuchi=0
        for j in range(len(start_probe_arrayf)):
            if (((start_probe_arrayf[j]>=start_array1[i]-10000 and end_probe_arrayf[j]<=start_array1[i]) or ( (start_probe_arrayf[j]<=start_array1[i]-10000 and end_probe_arrayf[j]>= start_array1[i]-10000) or (start_probe_arrayf[j]<=start_array1[i] and end_probe_arrayf[j]>=start_array1[i])))) and len(temp_array_ref0)<6: 
                temp_array_ref0.append(value_arrayf[j])
            if (((start_probe_arrayf[j]>=end_array1[i] and end_probe_arrayf[j]<=end_array1[i]+10000) or ( (start_probe_arrayf[j]<=end_array1[i] and end_probe_arrayf[j]>= end_array1[i]) or (start_probe_arrayf[j]<=end_array1[i]+10000 and end_probe_arrayf[j]>=end_array1[i]+10000)))) and len(temp_array_ref1)<6:
                temp_array_ref1.append(value_arrayf[j])
            
            if (((start_probe_arrayf[j]>=start_array1[i] and end_probe_arrayf[j]<=end_array1[i]) or ( (start_probe_arrayf[j]<=start_array1[i] and end_probe_arrayf[j]>= start_array1[i]) or (start_probe_arrayf[j]<=end_array1[i] and end_probe_arrayf[j]>=end_array1[i])))):
                #print start_array1[i]
                #print end_array1[i]
                #print start_probe_arrayf[j]
                #print end_probe_arrayf[j]
                if ( (start_probe_arrayf[j]<start_array1[i] and (end_probe_arrayf[j]-start_array1[i])>=(((end_probe_arrayf[j]-start_probe_arrayf[j])*30)/100)) or (end_probe_arrayf[j]>end_array1[i] and ((end_array1[i]-start_probe_arrayf[j])>=(((end_probe_arrayf[j]-start_probe_arrayf[j])*30)/100)))): 
                
                    countTf=countTf+1
                    temp_array_st.append(start_probe_arrayf[j])
                    temp_array_end.append(end_probe_arrayf[j])
                    temp_array_ref.append(value_arrayf[j])
                    #print 'entro'
                    if status_arrayf[j]=='C':
                        countCf=countCf+1
                    if status_arrayf[j]=='N':
                        countNf=countNf+1
                
                    controlla_nome=nome_probe_arrayf[j]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGenef=countGenef+1
                if (start_probe_arrayf[j]>=start_array1[i]) and end_probe_arrayf[j]<=end_array1[i]:
                    countTf=countTf+1
                    temp_array_ref.append(value_arrayf[j])
                    temp_array_st.append(start_probe_arrayf[j])
                    temp_array_end.append(end_probe_arrayf[j])
                
                    if status_arrayf[j]=='C':
                        countCf=countCf+1
                    if status_arrayf[j]=='N':
                        countNf=countNf+1
                
                    controlla_nome=nome_probe_arrayf[j]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGenef=countGenef+1
        
        for jj in range(len(start_array_gbk_f)):
            if (((start_array_gbk_f[jj]>=start_array1[i] and end_array_gbk_f[jj]<=end_array1[i]) or ( (start_array_gbk_f[jj]<=start_array1[i] and end_array_gbk_f[jj]>= start_array1[i]) or (start_array_gbk_f[jj]<=end_array1[i] and end_array_gbk_f[jj]>=end_array1[i]) or (start_array_gbk_f[jj]<=start_array1[i] and end_array_gbk_f[jj]>=end_array1[i]) ) )):
                
                #print start_array1[i]
                #print end_array1[i]
                #print start_probe_arrayf[j]
                #print end_probe_arrayf[j]
                if ( (start_array_gbk_f[jj]<start_array1[i] and (end_array_gbk_f[jj]-start_array1[i])>=(((end_array_gbk_f[jj]-start_array_gbk_f[jj])*10)/100)) or (end_array_gbk_f[jj]>end_array1[i] and ((end_array1[i]-start_array_gbk_f[jj])>=(((end_array_gbk_f[jj]-start_array_gbk_f[jj])*10)/100))) or (start_array_gbk_f[jj]>=start_array1[i] and end_array_gbk_f[jj]<=end_array1[i])): 
                    countGenef=countGenef+1
                    len_gbk=end_array_gbk_f[jj]-start_array_gbk_f[jj]
                    #len_sad=end_array1[i]-start_array1[i]
                    if start_array_gbk_f[jj]>=start_array1[i]:
                        if end_array_gbk_f[jj]>=end_array1[i]:
                            lenght=end_array1[i]-start_array_gbk_f[jj]
                        if end_array_gbk_f[jj]<=end_array1[i]:
                            lenght=end_array_gbk_f[jj]-start_array_gbk_f[jj]
                    if start_array_gbk_f[jj]<=start_array1[i]:
                        if end_array_gbk_f[jj]<=end_array1[i]:
                            
                            lenght=end_array_gbk_f[jj]-start_array1[i]
                        if end_array_gbk_f[jj]>=end_array1[i]:
                            lenght=end_array1[i]-start_array1[i]
                    perc=(len_gbk*30)/100
                    
                    if lenght>=perc:
                        count_gene_f=count_gene_f+1
                        #print count_gene_f
                        #print len(temp_array_end)
                        for hh in range(len(temp_array_end)):
                            #print hh
                            if temp_array_end[hh]>=start_array_gbk_f[jj]-60 and temp_array_st[hh]<=start_array_gbk_f[jj]-60 : 
                                count_nobuchi=count_nobuchi+1
                                #print count_nobuchi
                            if temp_array_st[hh]<=end_array_gbk_f[jj]+60 and temp_array_end[hh]>=end_array_gbk_f[jj]+60:
                                count_nobuchi=count_nobuchi+1
                                #print count_nobuchi
                    """
                    file_out.write("end_probe %s >= start_gbk %s;   start_probe %s <= end_gbk %s" %(temp_array_end[len(temp_array_end)-1], start_array_gbk_f[jj], temp_array_st[0], end_array_gbk_f[jj]))
                    if temp_array_end[len(temp_array_end)-1]>=start_array_gbk_f[jj] or temp_array_st[0]<=end_array_gbk_f[jj]:
                        count_nobuchi=count_nobuchi+1
                        file_out.write("nobuchi %s\n" %(count_nobuchi))
                    """
                        
        for z in range(len(start_probe_arrayr)):    
            if (((start_probe_arrayr[z]>=start_array1[i] and end_probe_arrayr[z]<=end_array1[i]) or ( (start_probe_arrayr[z]<=start_array1[i] and end_probe_arrayr[z]>= start_array1[i]) or (start_probe_arrayr[z]<=end_array1[i] and end_probe_arrayr[z]>=end_array1[i]))))  : 
                if ( (start_probe_arrayr[z]<start_array1[i] and (end_probe_arrayr[z]-start_array1[i])>=(((end_probe_arrayr[z]-start_probe_arrayr[z])*30)/100)) or (end_probe_arrayr[z]>end_array1[i] and ((end_array1[i]-start_probe_arrayr[z])>=(((end_probe_arrayr[z]-start_probe_arrayr[z])*30)/100)))): 
                
                    countTr=countTr+1
                    temp_array_refr.append(value_arrayr[z])
                    if status_arrayr[z]=='C':
                        countCr=countCr+1
                    if status_arrayr[z]=='N':
                        countNr=countNr+1
                    controlla_nome=nome_probe_arrayr[z]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGener=countGener+1
                if (start_probe_arrayr[z]>=start_array1[i]) and end_probe_arrayr[z]<=end_array1[i]:
                    countTr=countTr+1
                    temp_array_refr.append(value_arrayr[z])
                    if status_arrayr[z]=='C':
                        countCr=countCr+1
                    if status_arrayr[z]=='N':
                        countNr=countNr+1
                    controlla_nome=nome_probe_arrayr[z]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGener=countGener+1
                        
                        
        for zz in range(len(start_array_gbk_r)):
            if (((start_array_gbk_r[zz]>=start_array1[i] and end_array_gbk_r[zz]<=end_array1[i]) or ( (start_array_gbk_r[zz]<=start_array1[i] and end_array_gbk_r[zz]>= start_array1[i]) or (start_array_gbk_r[zz]<=end_array1[i] and end_array_gbk_r[zz]>=end_array1[i]) or (start_array_gbk_r[zz]<=start_array1[i] and end_array_gbk_r[zz]>=end_array1[i]) ))):
                
                #print start_array1[i]
                #print end_array1[i]
                #print start_probe_arrayf[j]
                #print end_probe_arrayf[j]
                if ( (start_array_gbk_r[zz]<start_array1[i] and (end_array_gbk_r[zz]-start_array1[i])>=(((end_array_gbk_r[zz]-start_array_gbk_r[zz])*10)/100)) or (end_array_gbk_r[zz]>end_array1[i] and ((end_array1[i]-start_array_gbk_r[zz])>=(((end_array_gbk_r[zz]-start_array_gbk_r[zz])*10)/100))) or (start_array_gbk_r[zz]>=start_array1[i] and end_array_gbk_r[zz]<=end_array1[i])): 
                    len_gbk=end_array_gbk_r[zz]-start_array_gbk_r[zz]
                    #len_sad=end_array1[i]-start_array1[i]
                    countGener=countGener+1
                    if start_array_gbk_r[zz]>=start_array1[i]:
                        if end_array_gbk_r[zz]>=end_array1[i]:
                            lenght=end_array1[i]-start_array_gbk_r[zz]
                        if end_array_gbk_r[zz]<=end_array1[i]:
                            lenght=end_array_gbk_r[zz]-start_array_gbk_r[zz]
                    if start_array_gbk_r[zz]<=start_array1[i]:
                        if end_array_gbk_r[zz]<=end_array1[i]:
                            lenght=end_array_gbk_r[zz]-start_array1[i]
                        if end_array_gbk_r[zz]>=end_array1[i]:                         
                            lenght=end_array1[i]-start_array1[i]
                    perc=(len_gbk*30)/100
                    if lenght>=perc:         
                        count_gene_r=count_gene_r+1
                    #print count_gene_r           
        
        
        
        
        #print temp_array_ref
        #print temp_array_ref0
        #print temp_array_ref1
        #print temp_array_refr   
        #file_utr.write("T %s, C %s, N %s\n" %(countTf, countCf, countNf))
        #if countCf+countNf==countTf and countCf>=2 and countNf>=3 and countCr==0:
        #if countCf+countNf==countTf and countCf>=2 and countNf>=3:
        
        #length_ref=len(temp_array_ref)
        #if len(temp_array_refr)==0:
         #   while k<length_ref:
          #      temp_array_refr.append(1.0)
           #     k=k+1
           
        #QUI   
         
        #if countCf>=2 and countNf>=3 and (max_array1[i]>=1 or max_array1[i]<=-1) and countCr==0:
        #if count_nobuchi!=0:
         #   print count_nobuchi
        
        if countCf>0 and count_gene_f==1 and countNf>=1 and  count_gene_r==0 and count_nobuchi>0 and (max_array1[i]>=thup or max_array1[i]<=thdown) and pvalue_array1[i]<=0.05: #
            #print 'ci sono'
            start_utr.append(start_array1[i])
            end_utr.append(end_array1[i])
            max_utr.append(max_array1[i])
            kutr=kutr+1
            strand='F'
            for cont in range(len(start_array_gbk_f)):
                if (((start_array_gbk_f[cont]>=start_array1[i] and end_array_gbk_f[cont]<=end_array1[i]) or ( (start_array_gbk_f[cont]<=start_array1[i] and end_array_gbk_f[cont]>= start_array1[i]) or (start_array_gbk_f[cont]<=end_array1[i] and end_array_gbk_f[cont]>=end_array1[i]) or (start_array_gbk_f[cont]<=start_array1[i] and end_array_gbk_f[cont]>=end_array1[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    file_utrf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, name_array_gbk_f[cont], start_array1[i], end_array1[i], end_array1[i]-start_array1[i], max_array1[i], ttest_array1[i], pvalue_array1[i]))
        
            #file_utrf.write("st %s  ---  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
       
        
        #if countCf>=2 and countNf>=3 and (max_array1[i]>=1 or max_array1[i]<=-1) and countCr>0:
        if countCf>0 and count_gene_f==1 and countNf>=1  and count_gene_r>0 and  count_nobuchi>0  and (max_array1[i]>=thup or max_array1[i]<=thdown) and pvalue_array1[i]<=0.05: #
            start_utr_over.append(start_array1[i])
            end_utr_over.append(end_array1[i])
            max_utr_over.append(max_array1[i])
            kutrover=kutrover+1
            strand='F'
            for cont in range(len(start_array_gbk_f)):
                if (((start_array_gbk_f[cont]>=start_array1[i] and end_array_gbk_f[cont]<=end_array1[i]) or ( (start_array_gbk_f[cont]<=start_array1[i] and end_array_gbk_f[cont]>= start_array1[i]) or (start_array_gbk_f[cont]<=end_array1[i] and end_array_gbk_f[cont]>=end_array1[i]) or (start_array_gbk_f[cont]<=start_array1[i] and end_array_gbk_f[cont]>=end_array1[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    gene=name_array_gbk_f[cont]
            for contr in range(len(start_array_gbk_r)):
                if (((start_array_gbk_r[contr]>=start_array1[i] and end_array_gbk_r[contr]<=end_array1[i]) or ( (start_array_gbk_r[contr]<=start_array1[i] and end_array_gbk_r[contr]>= start_array1[i]) or (start_array_gbk_r[contr]<=end_array1[i] and end_array_gbk_r[contr]>=end_array1[i]) or (start_array_gbk_r[contr]<=start_array1[i] and end_array_gbk_r[contr]>=end_array1[i]) ))):              
                    orfover=orfover+name_array_gbk_r[contr]
            file_utr_overf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, gene, start_array1[i], end_array1[i], end_array1[i]-start_array1[i], orfover , max_array1[i], ttest_array1[i], pvalue_array1[i]))
        
            orfover=''
            #file_utr_overf.write("st %s  ---  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
        
        
        
        #if (countNf==countTf and (max_array1[i]>=1 or max_array1[i]<=-1) and countCr!=0 and countNr<3):
        if (count_gene_f==0  and count_gene_r>0) and (max_array1[i]>=thup or max_array1[i]<=thdown) and pvalue_array1[i]<=0.05:
            start_anti.append(start_array1[i])
            end_anti.append(end_array1[i])
            max_anti.append(max_array1[i])
            kanti=kanti+1
            strand='F'
            for cont in range(len(start_array_gbk_r)):
                if (((start_array_gbk_r[cont]>=start_array1[i] and end_array_gbk_r[cont]<=end_array1[i]) or ( (start_array_gbk_r[cont]<=start_array1[i] and end_array_gbk_r[cont]>= start_array1[i]) or (start_array_gbk_r[cont]<=end_array1[i] and end_array_gbk_r[cont]>=end_array1[i]) or (start_array_gbk_r[cont]<=start_array1[i] and end_array_gbk_r[cont]>=end_array1[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    file_antisensef.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, name_array_gbk_r[cont], start_array1[i], end_array1[i], end_array1[i]-start_array1[i], max_array1[i], ttest_array1[i], pvalue_array1[i]))
        
        
        
        
        #if countGenef==((2*countTf)/3) and countCf>=6 and (max_array1[i]>=1 or max_array1[i]<=-1):
        if  count_gene_f>=2 and (max_array1[i]>=thup or max_array1[i]<=thdown) and pvalue_array1[i]<=0.05:
            start_operonsf.append(start_array1[i])
            end_operonsf.append(end_array1[i])
            max_operonsf.append(max_array1[i])
            kop=kop+1
            strand='F'
            for cont in range(len(start_array_gbk_f)):
                if (((start_array_gbk_f[cont]>=start_array1[i] and end_array_gbk_f[cont]<=end_array1[i]) or ( (start_array_gbk_f[cont]<=start_array1[i] and end_array_gbk_f[cont]>= start_array1[i]) or (start_array_gbk_f[cont]<=end_array1[i] and end_array_gbk_f[cont]>=end_array1[i]) or (start_array_gbk_f[cont]<=start_array1[i] and end_array_gbk_f[cont]>=end_array1[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    orfs=orfs+name_array_gbk_f[cont]
            file_operonsf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, orfs, start_array1[i], end_array1[i], end_array1[i]-start_array1[i], max_array1[i], ttest_array1[i], pvalue_array1[i]))
        
            orfs=''
        
            #file_operonsf.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
        
         #QUI
       
        #print thup
        #print thdown
        #if countCf==0 and countCr==0 and (max_array1[i]>=1 or max_array1[i]<=-1) and countTr>0   and (end_array1[i]-start_array1[i])<800: #and stats.ttest_ind(temp_array_ref, temp_array_refr)[1]<0.05: #and  stats.ttest_ind(max_array1[i], max_array1[i+1], 0)[1]<0.05: #and abs(max_array1[i])/2>abs(max_array1[i+1]) and abs(max_array1[i])/2>abs(max_array1[i-1]) and abs(max_array1[i])/2>abs(max_array2[i]):
        if countGenef==0 and countGener==0 and  (end_array1[i]-start_array1[i])<800 and (max_array1[i]>=thup or max_array1[i]<=thdown):# and pvalue_array1[i]<=0.05:    #and countTr>0:  # #and countTr>0
            start_interf.append(start_array1[i])
            end_interf.append(end_array1[i])
            max_interf.append(max_array1[i])
            strand='F'
            #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
            #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
            #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
            kinter=kinter+1
            for cont in range(len(start_array_gbk_glob)):
                if start_array1[i]>=start_array_gbk_glob[cont]:
                    temp_glob_st.append(name_array_gbk_glob[cont])
                    temp_glob_strand_st.append(strand_glob_array[cont])
                if end_array1[i]<=end_array_gbk_glob[cont]:
                    temp_glob_end.append(name_array_gbk_glob[cont])
                    temp_glob_strand_end.append(strand_glob_array[cont])
            if len(temp_glob_st)!=0 and len(temp_glob_end)!=0 and temp_glob_st[len(temp_glob_st)-1]!=temp_glob_end[0]:
                prima=temp_glob_st[len(temp_glob_st)-1]
                strandprima=temp_glob_strand_st[len(temp_glob_strand_st)-1]
                dopo=temp_glob_end[0]
                stranddopo=temp_glob_strand_end[0]
                file_interf.write("%s\t%s\t%s\t%s>%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, prima, dopo, strandprima, stranddopo, start_array1[i], end_array1[i], end_array1[i]-start_array1[i], max_array1[i], ttest_array1[i], pvalue_array1[i]))
        
            temp_glob_st=[]
            temp_glob_end=[]
            temp_glob_strand_st=[]
            temp_glob_strand_end=[]
                
                #if (((start_array_gbk_glob[cont]>=start_array1[i] and end_array_gbk_glob[cont]<=end_array1[i]) or ( (start_array_gbk_glob[cont]<=start_array1[i] and end_array_gbk_glob[cont]>= start_array1[i]) or (start_array_gbk_glob[cont]<=end_array1[i] and end_array_gbk_glob[cont]>=end_array1[i]) or (start_array_gbk_glob[cont]<=start_array1[i] and end_array_gbk_glob[cont]>=end_array1[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    #file_interf.write("%s\t%s\t%s\t%s\t%s\t\n" %(strand, name_array_gbk_glob[cont], start_array1[i], end_array1[i], end_array1[i]-start_array1[i]))
            #file_interf.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
            
            #file_interf.write("prima st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
            #file_interf.write("len %s\n " %(len(temp_array_ref)))
            #print temp_array_ref
            #print temp_array_ref0
            #print temp_array_ref1
            #print temp_array_refr
            #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
            #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
            #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
            #print temp_array_ref
            #print temp_array_refr
            """
            if len(temp_array_ref)>1:
                if (stats.ttest_ind(temp_array_ref, temp_array_refr)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref0)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref1)[1]<0.05) or (max_array1[i]>=3 or max_array1[i]<=-3) :       
                    start_interf.append(start_array1[i])
                    end_interf.append(end_array1[i])
                    max_interf.append(max_array1[i])
                    #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                    #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                    #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
                    kinter=kinter+1
                    file_interf.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
            if len(temp_array_ref)==1:
                if (abs(temp_array_ref[0])/2>abs(temp_array_refr[0]) and stats.ttest_ind(temp_array_ref, temp_array_ref0)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref1)[1]<0.05) or (max_array1[i]>=3 or max_array1[i]<=-3):       
                    start_interf.append(start_array1[i])
                    end_interf.append(end_array1[i])
                    max_interf.append(max_array1[i])
                    #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                    #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                    #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
                    kinter=kinter+1
                    file_interf.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
        
            """      
                    
        temp_array_ref=[]
        temp_array_refr=[]
        temp_array_ref0=[]
        temp_array_ref1=[]
        temp_array_st=[]
        temp_array_end=[]
            
   
    file_out.write("F UTR %s\n" %(kutr))    
    file_out.write("F UTR over %s\n" %(kutrover)) 
    file_out.write("F antisense %s\n" %(kanti))
    file_out.write("F operons %s\n" %(kop))
    file_out.write("F intergenic %s\n" %(kinter))   
    print "F strand done"
    
    
    i=0
    j=0
    z=0
    jj=0
    hh=0
    zz=0
    count_gene_r=0
    count_gene_f=0
    start_utrr=[]
    end_utrr=[]
    max_utrr=[]
    start_utr_overr=[]
    end_utr_overr=[]
    max_utr_overr=[]
    start_antir=[]
    end_antir=[]
    max_antir=[]
    start_operonsr=[]
    end_operonsr=[]
    max_operonsr=[]
    start_interr=[]
    end_interr=[]
    max_interr=[]
    temp_array_ref=[]
    temp_array_refr=[]
    temp_array_ref1=[]
    temp_array_ref0=[]
    temp_array_st=[]
    temp_array_end=[]
    temp_glob_st=[]
    temp_glob_end=[]
    temp_glob_strand_st=[]
    temp_glob_strand_end=[]
    temp_array_st=[]
    temp_array_end=[]
    orfs=''
    orfover=''
    kutr=0
    kutrover=0
    kanti=0
    kop=0
    kinter=0
    count_nobuchir=0
    for i in range(len(start_array2)):
        countCf=0
        countCr=0
        countNf=0
        countNr=0
        countTf=0
        countTr=0
        
        count_nobuchir=0
        count_gene_r=0
        count_gene_f=0
        
        countGenef=0
        countGener=0
        for j in range(len(start_probe_arrayr)):
            if (((start_probe_arrayr[j]>=start_array2[i]-10000 and end_probe_arrayr[j]<=start_array2[i]) or ( (start_probe_arrayr[j]<=start_array2[i]-10000 and end_probe_arrayr[j]>= start_array2[i]-10000) or (start_probe_arrayr[j]<=start_array2[i] and end_probe_arrayr[j]>=start_array2[i])))) and len(temp_array_ref0)<6: 
                temp_array_ref0.append(value_arrayr[j])
            if (((start_probe_arrayr[j]>=end_array2[i] and end_probe_arrayr[j]<=end_array2[i]+10000) or ( (start_probe_arrayr[j]<=end_array2[i] and end_probe_arrayr[j]>= end_array2[i]) or (start_probe_arrayr[j]<=end_array2[i]+10000 and end_probe_arrayr[j]>=end_array2[i]+10000)))) and len(temp_array_ref1)<6:
                temp_array_ref1.append(value_arrayr[j])
            
            if (((start_probe_arrayr[j]>=start_array2[i] and end_probe_arrayr[j]<=end_array2[i]) or ( (start_probe_arrayr[j]<=start_array2[i] and end_probe_arrayr[j]>= start_array2[i]) or (start_probe_arrayr[j]<=end_array2[i] and end_probe_arrayr[j]>=end_array2[i])))):
                #print start_array1[i]
                #print end_array1[i]
                #print start_probe_arrayf[j]
                #print end_probe_arrayf[j]
                if ( (start_probe_arrayr[j]<start_array2[i] and (end_probe_arrayr[j]-start_array2[i])>=(((end_probe_arrayr[j]-start_probe_arrayr[j])*30)/100)) or (end_probe_arrayr[j]>end_array2[i] and ((end_array2[i]-start_probe_arrayr[j])>=(((end_probe_arrayr[j]-start_probe_arrayr[j])*30)/100)))): 
                
                    countTr=countTr+1
                    temp_array_ref.append(value_arrayr[j])
                    temp_array_st.append(start_probe_arrayr[j])
                    temp_array_end.append(end_probe_arrayr[j])
                
                    if status_arrayr[j]=='C':
                        countCr=countCr+1
                    if status_arrayr[j]=='N':
                        countNr=countNr+1
                
                    controlla_nome=nome_probe_arrayr[j]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGener=countGener+1
                if (start_probe_arrayr[j]>=start_array2[i]) and end_probe_arrayr[j]<=end_array2[i]:
                    countTr=countTr+1
                    temp_array_ref.append(value_arrayr[j])
                    temp_array_st.append(start_probe_arrayr[j])
                    temp_array_end.append(end_probe_arrayr[j])
                
                    if status_arrayr[j]=='C':
                        countCr=countCr+1
                    if status_arrayr[j]=='N':
                        countNr=countNr+1
                
                    controlla_nome=nome_probe_arrayr[j]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGener=countGener+1
        
        
        for jj in range(len(start_array_gbk_r)):
            if (((start_array_gbk_r[jj]>=start_array2[i] and end_array_gbk_r[jj]<=end_array2[i]) or ( (start_array_gbk_r[jj]<=start_array2[i] and end_array_gbk_r[jj]>= start_array2[i]) or (start_array_gbk_r[jj]<=end_array2[i] and end_array_gbk_r[jj]>=end_array2[i]) or (start_array_gbk_r[jj]<=start_array2[i] and end_array_gbk_r[jj]>=end_array2[i]) ) )):
                
                #print start_array1[i]
                #print end_array1[i]
                #print start_probe_arrayf[j]
                #print end_probe_arrayf[j]
                if ( (start_array_gbk_r[jj]<start_array2[i] and (end_array_gbk_r[jj]-start_array2[i])>=(((end_array_gbk_r[jj]-start_array_gbk_r[jj])*10)/100)) or (end_array_gbk_r[jj]>end_array2[i] and ((end_array2[i]-start_array_gbk_r[jj])>=(((end_array_gbk_r[jj]-start_array_gbk_r[jj])*10)/100))) or (start_array_gbk_r[jj]>=start_array2[i] and end_array_gbk_r[jj]<=end_array2[i])): 
                    len_gbk=end_array_gbk_r[jj]-start_array_gbk_r[jj]
                    #len_sad=end_array2[i]-start_array2[i]
                    
                    
                    if start_array_gbk_r[jj]>=start_array2[i]:
                        if end_array_gbk_r[jj]>=end_array2[i]:
                            lenght=end_array2[i]-start_array_gbk_r[jj]
                        if end_array_gbk_r[jj]<=end_array2[i]:
                            lenght=end_array_gbk_r[jj]-start_array_gbk_r[jj]
                    if start_array_gbk_r[jj]<=start_array2[i]:
                        if end_array_gbk_r[jj]<=end_array2[i]:
                            lenght=end_array_gbk_r[jj]-start_array2[i]
                        if end_array_gbk_r[jj]>=end_array2[i]:
                            lenght=end_array2[i]-start_array2[i]
                    
                    
                    perc=(len_gbk*30)/100
                    countGener=countGener+1
                    if lenght>=perc:
                        count_gene_r=count_gene_r+1
                        for hh in range(len(temp_array_end)):
                            if temp_array_end[hh]>=start_array_gbk_r[jj]-60 and temp_array_st[hh]<=start_array_gbk_r[jj]-60 : 
                                count_nobuchir=count_nobuchir+1
                            if temp_array_st[hh]<=end_array_gbk_r[jj]+60 and temp_array_end[hh]>=end_array_gbk_r[jj]+60:
                                count_nobuchir=count_nobuchir+1
        
        
                       
        for z in range(len(start_probe_arrayf)):
            
            if (((start_probe_arrayf[z]>=start_array2[i] and end_probe_arrayf[z]<=end_array2[i]) or ( (start_probe_arrayf[z]<=start_array2[i] and end_probe_arrayf[z]>= start_array2[i]) or (start_probe_arrayf[z]<=end_array2[i] and end_probe_arrayf[z]>=end_array2[i]))))  : 
                
                if ( (start_probe_arrayf[z]<start_array2[i] and (end_probe_arrayf[z]-start_array2[i])>=(((end_probe_arrayf[z]-start_probe_arrayf[z])*30)/100)) or (end_probe_arrayf[z]>end_array2[i] and ((end_array2[i]-start_probe_arrayf[z])>=(((end_probe_arrayf[z]-start_probe_arrayf[z])*30)/100)))): 
                    
                    countTf=countTf+1
                    temp_array_refr.append(value_arrayf[z])
                    if status_arrayf[z]=='C':
                        countCf=countCf+1
                    if status_arrayf[z]=='N':
                        countNf=countNf+1
                    controlla_nome=nome_probe_arrayf[z]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGenef=countGenef+1
                if (start_probe_arrayf[z]>=start_array2[i]) and end_probe_arrayf[z]<=end_array2[i]:
                    countTf=countTf+1
                    temp_array_refr.append(value_arrayf[z])
                    if status_arrayf[z]=='C':
                        countCf=countCf+1
                    if status_arrayf[z]=='N':
                        countNf=countNf+1
                    controlla_nome=nome_probe_arrayf[z]
                    if fnmatch.fnmatch(controlla_nome, '*_B'):
                        countGenef=countGenef+1
        
        
        
        for zz in range(len(start_array_gbk_f)):
            if (((start_array_gbk_f[zz]>=start_array2[i] and end_array_gbk_f[zz]<=end_array2[i]) or ( (start_array_gbk_f[zz]<=start_array2[i] and end_array_gbk_f[zz]>= start_array2[i]) or (start_array_gbk_f[zz]<=end_array2[i] and end_array_gbk_f[zz]>=end_array2[i]) or (start_array_gbk_f[zz]<=start_array2[i] and end_array_gbk_f[zz]>=end_array2[i]) ))):
                
                ##print start_array1[i]
                #print end_array1[i]
                #print start_probe_arrayf[j]
                #print end_probe_arrayf[j]
                if ( (start_array_gbk_f[zz]<start_array2[i] and (end_array_gbk_f[zz]-start_array2[i])>=(((end_array_gbk_f[zz]-start_array_gbk_f[zz])*10)/100)) or (end_array_gbk_f[zz]>end_array2[i] and ((end_array2[i]-start_array_gbk_f[zz])>=(((end_array_gbk_f[zz]-start_array_gbk_f[zz])*10)/100))) or (start_array_gbk_f[zz]>=start_array2[i] and end_array_gbk_f[zz]<=end_array2[i])): 
                    len_gbk=end_array_gbk_f[zz]-start_array_gbk_f[zz]
                    #len_sad=end_array2[i]-start_array2[i]
                    
                    
                    if start_array_gbk_f[zz]>=start_array2[i]:
                        if end_array_gbk_f[zz]>=end_array2[i]:
                            lenght=end_array2[i]-start_array_gbk_f[zz]
                        if end_array_gbk_f[zz]<=end_array2[i]:
                            lenght=end_array_gbk_f[zz]-start_array_gbk_f[zz]
                    if start_array_gbk_f[zz]<=start_array2[i]:
                        if end_array_gbk_f[zz]<=end_array2[i]:
                            
                            lenght=end_array_gbk_f[zz]-start_array2[i]
                        if end_array_gbk_f[zz]>=end_array2[i]:
                            lenght=end_array2[i]-start_array2[i]
                            
                    
                    
                    perc=(len_gbk*30)/100
                    countGenef=countGenef+1
                    if lenght>=perc:
                        count_gene_f=count_gene_f+1
                    
                    #print count_gene_f          
        
        
        
            
        #file_utr.write("T %s, C %s, N %s\n" %(countTf, countCf, countNf))
        #if countCf+countNf==countTf and countCf>=2 and countNf>=3 and countCr==0:
        #if countCf+countNf==countTf and countCf>=2 and countNf>=3:
       
         #QUI
        #if countCr>=2 and countNr>=3 and (max_array2[i]>=1 or max_array2[i]<=-1) and countCf==0:
        #if countCf>0 and count_gene_f==1 and countNf>=1 and (max_array1[i]>=1 or max_array1[i]<=-1) and count_gene_r==0:
        if countCr>0 and count_gene_r==1 and countNr>=1 and (max_array2[i]>=thup or max_array2[i]<=thdown) and count_gene_f==0 and pvalue_array2[i]<=0.05 and count_nobuchir>0:
            start_utrr.append(start_array2[i])
            end_utrr.append(end_array2[i])
            max_utrr.append(max_array2[i])
            kutr=kutr+1
            #file_utrr.write("st %s  ---  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
            strand='R'
            for cont in range(len(start_array_gbk_r)):
                if (((start_array_gbk_r[cont]>=start_array2[i] and end_array_gbk_r[cont]<=end_array2[i]) or ( (start_array_gbk_r[cont]<=start_array2[i] and end_array_gbk_r[cont]>= start_array2[i]) or (start_array_gbk_r[cont]<=end_array2[i] and end_array_gbk_r[cont]>=end_array2[i]) or (start_array_gbk_r[cont]<=start_array2[i] and end_array_gbk_r[cont]>=end_array2[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    file_utrr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, name_array_gbk_r[cont], start_array2[i], end_array2[i], end_array2[i]-start_array2[i], max_array2[i], ttest_array2[i], pvalue_array2[i]))
        
        
        
            
        
        #if countCr>=2 and countNr>=3 and (max_array2[i]>=1 or max_array2[i]<=-1) and countCf>0:
        #if countCf>0 and count_gene_f==1 and countNf>=1 and (max_array1[i]>=1 or max_array1[i]<=-1) and count_gene_r>0:
        if countCr>0 and count_gene_r==1 and countNr>=1 and count_gene_f>0 and count_nobuchir>0 and (max_array2[i]>=thup or max_array2[i]<=thdown) and pvalue_array2[i]<=0.05: #
            start_utr_overr.append(start_array2[i])
            end_utr_overr.append(end_array2[i])
            max_utr_overr.append(max_array2[i])
            kutrover=kutrover+1
            #file_utr_overr.write("st %s  ---  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
            strand='R'
            for cont in range(len(start_array_gbk_r)):
                if (((start_array_gbk_r[cont]>=start_array2[i] and end_array_gbk_r[cont]<=end_array2[i]) or ( (start_array_gbk_r[cont]<=start_array2[i] and end_array_gbk_r[cont]>= start_array2[i]) or (start_array_gbk_r[cont]<=end_array2[i] and end_array_gbk_r[cont]>=end_array2[i]) or (start_array_gbk_r[cont]<=start_array2[i] and end_array_gbk_r[cont]>=end_array2[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    gene=name_array_gbk_r[cont]
            for contr in range(len(start_array_gbk_f)):
                if (((start_array_gbk_f[contr]>=start_array2[i] and end_array_gbk_f[contr]<=end_array2[i]) or ( (start_array_gbk_f[contr]<=start_array2[i] and end_array_gbk_f[contr]>= start_array2[i]) or (start_array_gbk_f[contr]<=end_array2[i] and end_array_gbk_f[contr]>=end_array2[i]) or (start_array_gbk_f[contr]<=start_array2[i] and end_array_gbk_f[contr]>=end_array2[i]) ))):              
                    orfover=orfover+name_array_gbk_f[contr]
            file_utr_overr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, gene, start_array2[i], end_array2[i], end_array2[i]-start_array2[i], orfover, max_array2[i], ttest_array2[i], pvalue_array2[i]))
        
            orfover=''
        
        
        
        #if (countNr==countTr and (max_array2[i]>=1 or max_array2[i]<=-1) and countCf!=0 and countNf<3):
        #if (count_gene_f==0 and (max_array1[i]>=1 or max_array1[i]<=-1) and count_gene_r>0):
        if (count_gene_r==0 and  count_gene_f>0) and (max_array2[i]>=thup or max_array2[i]<=thdown) and  pvalue_array2[i]<=0.05: 
            start_antir.append(start_array2[i])
            end_antir.append(end_array2[i])
            max_antir.append(max_array2[i])
            kanti=kanti+1
            #file_antisenser.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
            strand='R'
            for cont in range(len(start_array_gbk_f)):
                if (((start_array_gbk_f[cont]>=start_array2[i] and end_array_gbk_f[cont]<=end_array2[i]) or ( (start_array_gbk_f[cont]<=start_array2[i] and end_array_gbk_f[cont]>= start_array2[i]) or (start_array_gbk_f[cont]<=end_array2[i] and end_array_gbk_f[cont]>=end_array2[i]) or (start_array_gbk_f[cont]<=start_array2[i] and end_array_gbk_f[cont]>=end_array2[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    file_antisenser.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, name_array_gbk_f[cont], start_array2[i], end_array2[i], end_array2[i]-start_array2[i], max_array2[i], ttest_array2[i], pvalue_array2[i]))
        
        
        
        
        
        
        #if countGener==((2*countTr)/3) and countCr>=6 and (max_array2[i]>=1 or max_array2[i]<=-1):
        #if  count_gene_f>=2 and (max_array1[i]>=1 or max_array1[i]<=-1):
        if  count_gene_r>=2 and (max_array2[i]>=thup or max_array2[i]<=thdown) and pvalue_array2[i]<=0.05:
            start_operonsr.append(start_array2[i])
            end_operonsr.append(end_array2[i])
            max_operonsr.append(max_array2[i])
            kop=kop+1
            #file_operonsr.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
            strand='R'
            for cont in range(len(start_array_gbk_r)):
                if (((start_array_gbk_r[cont]>=start_array2[i] and end_array_gbk_r[cont]<=end_array2[i]) or ( (start_array_gbk_r[cont]<=start_array2[i] and end_array_gbk_r[cont]>= start_array2[i]) or (start_array_gbk_r[cont]<=end_array2[i] and end_array_gbk_r[cont]>=end_array2[i]) or (start_array_gbk_r[cont]<=start_array2[i] and end_array_gbk_r[cont]>=end_array2[i]) ))):
                    #file_antisensef.write("st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
                    orfs=orfs+name_array_gbk_r[cont]
            file_operonsr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, orfs, start_array2[i], end_array2[i], end_array2[i]-start_array2[i], max_array2[i], ttest_array2[i], pvalue_array2[i]))
        
            orfs=''
        
        
        #if len(temp_array_refr)==0:
         #   while len(temp_array_refr)<=len(temp_array_ref):
          #      temp_array_refr.append(0)
            
        
        #print temp_array_ref
        #print temp_array_ref0
        #print temp_array_ref1
        #print temp_array_refr
        
        
        #QUI
       
        #if countCf==0 and countCr==0 and (max_array2[i]>=1 or max_array2[i]<=-1) and countTr>0   and (end_array2[i]-start_array2[i])<800: #and stats.ttest_ind(temp_array_ref, temp_array_refr)[1]<0.05: #and  stats.ttest_ind(max_array1[i], max_array1[i+1], 0)[1]<0.05: #and abs(max_array1[i])/2>abs(max_array1[i+1]) and abs(max_array1[i])/2>abs(max_array1[i-1]) and abs(max_array1[i])/2>abs(max_array2[i]):
        #if count_gene_f==0 and count_gene_r==0 and (max_array1[i]>=1 or max_array1[i]<=-1) and countTr>0   and (end_array1[i]-start_array1[i])<800:
        if countGener==0 and countGenef==0 and (end_array2[i]-start_array2[i])<800 and (max_array2[i]>=thup or max_array2[i]<=thdown) :#    and pvalue_array2[i]<=0.05: #and countTf>0
            strand='R'
            
            #file_interf.write("prima st %s  --  end %s  ---  valmax %s\n" %(start_array1[i], end_array1[i], max_array1[i]))
            #file_interf.write("len %s\n " %(len(temp_array_ref)))
            start_interr.append(start_array2[i])
            end_interr.append(end_array2[i])
            max_interr.append(max_array2[i])
            #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
            #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
            #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
            kinter=kinter+1
            #file_interr.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
            for cont in range(len(start_array_gbk_glob)):
                if start_array2[i]>=start_array_gbk_glob[cont]:
                    temp_glob_st.append(name_array_gbk_glob[cont])
                    temp_glob_strand_st.append(strand_glob_array[cont])
                if end_array2[i]<=end_array_gbk_glob[cont]:
                    temp_glob_end.append(name_array_gbk_glob[cont])
                    temp_glob_strand_end.append(strand_glob_array[cont])
            if len(temp_glob_st)!=0 and len(temp_glob_end)!=0 and temp_glob_st[len(temp_glob_st)-1]!=temp_glob_end[0]:
                prima=temp_glob_st[len(temp_glob_st)-1]
                strandprima=temp_glob_strand_st[len(temp_glob_strand_st)-1]
                dopo=temp_glob_end[0]
                stranddopo=temp_glob_strand_end[0]
                file_interr.write("%s\t%s\t%s\t%s<%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n" %(strand, prima, dopo, strandprima, stranddopo, start_array2[i], end_array2[i], end_array2[i]-start_array2[i], max_array2[i], ttest_array2[i], pvalue_array2[i]))
                
            temp_glob_st=[]
            temp_glob_end=[]
            temp_glob_strand_st=[]
            temp_glob_strand_end=[]
            
            """
            if len(temp_array_ref)>1:
                if len(temp_array_refr)!=0:
                    if (stats.ttest_ind(temp_array_ref, temp_array_refr)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref0)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref1)[1]<0.05) or  (max_array2[i]>=3 or max_array2[i]<=-3) :       
                        start_interr.append(start_array2[i])
                        end_interr.append(end_array2[i])
                        max_interr.append(max_array2[i])
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                        #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
                        kinter=kinter+1
                        file_interr.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
                else:
                    #print temp_array_ref
                    #print temp_array_ref0
                    #print temp_array_ref1
                    #print temp_array_refr
                    #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                    #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                    if (stats.ttest_ind(temp_array_ref, temp_array_ref0)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref1)[1]<0.05) or   (max_array2[i]>=3 or max_array2[i]<=-3):
                        #print "entro"
                        start_interr.append(start_array2[i])
                        end_interr.append(end_array2[i])
                        max_interr.append(max_array2[i])
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                        #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
                        kinter=kinter+1
                        file_interr.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
            else:
                if len(temp_array_refr)!=0:
                    if (abs(temp_array_ref[0])/2>abs(temp_array_refr[0]) and stats.ttest_ind(temp_array_ref, temp_array_ref0)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref1)[1]<0.05) or  (max_array2[i]>=3 or max_array2[i]<=-3):       
                        start_interr.append(start_array2[i])
                        end_interr.append(end_array2[i])
                        max_interr.append(max_array2[i])
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                        #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
                        kinter=kinter+1
                        file_interr.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
                else:
                    if (stats.ttest_ind(temp_array_ref, temp_array_ref0)[1]<0.05 and stats.ttest_ind(temp_array_ref, temp_array_ref1)[1]<0.05) or  (max_array2[i]>=3 or max_array2[i]<=-3):       
                        start_interr.append(start_array2[i])
                        end_interr.append(end_array2[i])
                        max_interr.append(max_array2[i])
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref0)
                        #print stats.ttest_ind(temp_array_ref, temp_array_ref1)
                        #print stats.ttest_ind(temp_array_ref, temp_array_refr)[1]
                        kinter=kinter+1
                        file_interr.write("st %s  --  end %s  ---  valmax %s\n" %(start_array2[i], end_array2[i], max_array2[i]))
                    
            """    
        temp_array_ref=[]
        temp_array_refr=[]
        temp_array_ref0=[]
        temp_array_ref1=[]
        temp_array_st=[]
        temp_array_end=[]
    
   
    file_out.write("R UTR %s\n" %(kutr))    
    file_out.write("R UTR over %s\n" %(kutrover)) 
    file_out.write("R antisense %s\n" %(kanti))
    file_out.write("R operons %s\n" %(kop))

    file_out.write("R intergenic %s\n" %(kinter))   
    print "R strand done"
    
            
    matricenew=[]
    pippo=0
    for pippo in range(len(start_utr)):
        array_internon=[]       
        array_internon0=start_utr[pippo]
        array_internon1=end_utr[pippo]
        array_internon2= max_utr[pippo]
        array_internon.append(array_internon0)
        array_internon.append(array_internon1)
        array_internon.append(array_internon2)
        matricenew.append(array_internon)
        
    #print matrice
    
    matricenew.sort()
    #matricenew.reverse()
    start_ord_utr=[]
    end_ord_utr=[]
    maxval_ord_utr=[]
    i=0
    for i in range(len(matricenew)):      
        start_ord_utr.append(matricenew[i][0])
        end_ord_utr.append(matricenew[i][1])
        maxval_ord_utr.append(matricenew[i][2])
    
    print "artemis files preparation F strand"
    if len(maxval_ord_utr)>0:
        max_av=max(maxval_ord_utr)
        min_av=min(maxval_ord_utr)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---UTR artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_utrf.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_utrf.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        for me in range(len(start_ord_utr)):
            medav=maxval_ord_utr[me]
            st=start_ord_utr[me]
            ed=end_ord_utr[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_utrf.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_utrf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_utrf.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_utrf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_utrf.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_utrf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_utrf.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_utrf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_utrf.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_utrf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_utrf.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_utrf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No UTR"
    
    
    
          
    matrice=[]
    pippo=0
    for pippo in range(len(start_utr_over)):
        array_interno=[]        
        array_interno0=start_utr_over[pippo]
        array_interno1=end_utr_over[pippo]
        array_interno2= max_utr_over[pippo]
        array_interno.append(array_interno0)
        array_interno.append(array_interno1)
        array_interno.append(array_interno2)
        matrice.append(array_interno)
            
            
    #print matrice
            
    matrice.sort()
    #matrice.reverse()
    start_ord_utr_over=[]
    end_ord_utr_over=[]
    maxval_ord_utr_over=[]
    i=0
    for i in range(len(matrice)):      
        start_ord_utr_over.append(matrice[i][0])
        end_ord_utr_over.append(matrice[i][1])
        maxval_ord_utr_over.append(matrice[i][2])
    
    
    
    if len(maxval_ord_utr_over)>0:
        max_av=max(maxval_ord_utr_over)
        min_av=min(maxval_ord_utr_over)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---overlapping UTR artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_utr_overf.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_utr_overf.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_utr_over)):
            medav=maxval_ord_utr_over[me]
            st=start_ord_utr_over[me]
            ed=end_ord_utr_over[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_utr_overf.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_utr_overf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_utr_overf.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_utr_overf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_utr_overf.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_utr_overf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_utr_overf.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_utr_overf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_utr_overf.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_utr_overf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_utr_overf.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_utr_overf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        #for zz in range(len(start_utr_over)):
        #    file_utr_over.write("st %s  ---  end %s\n" %(start_utr_over[zz], end_utr_over[zz]))
    
    else:
        print "No overlapping UTR"
    
    matricenewn=[]
    pippo=0
    for pippo in range(len(start_anti)):
        array_internonn=[]       
        array_internonn0=start_anti[pippo]
        array_internonn1=end_anti[pippo]
        array_internonn2= max_anti[pippo]
        array_internonn.append(array_internonn0)
        array_internonn.append(array_internonn1)
        array_internonn.append(array_internonn2)
        matricenewn.append(array_internonn)
        
        #print matrice
        
    matricenewn.sort()
    #matricenew.reverse()
    start_ord_anti=[]
    end_ord_anti=[]
    maxval_ord_anti=[]
    i=0
    for i in range(len(matricenewn)):      
        start_ord_anti.append(matricenewn[i][0])
        end_ord_anti.append(matricenewn[i][1])
        maxval_ord_anti.append(matricenewn[i][2])
    
    
    if len(maxval_ord_anti)>0:
        max_av=max(maxval_ord_anti)
        min_av=min(maxval_ord_anti)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Antisense artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_antisensef.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_antisensef.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_anti)):
            medav=maxval_ord_anti[me]
            st=start_ord_anti[me]
            ed=end_ord_anti[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_antisensef.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_antisensef.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_antisensef.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_antisensef.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_antisensef.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_antisensef.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_antisensef.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_antisensef.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_antisensef.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_antisensef.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_antisensef.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_antisensef.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No antisense"
    
    
    
    
    
    matricenewnop=[]
    pippo=0
    for pippo in range(len(start_operonsf)):
        array_internonnop=[]       
        array_internonnop0=start_operonsf[pippo]
        array_internonnop1=end_operonsf[pippo]
        array_internonnop2= max_operonsf[pippo]
        array_internonnop.append(array_internonnop0)
        array_internonnop.append(array_internonnop1)
        array_internonnop.append(array_internonnop2)
        matricenewnop.append(array_internonnop)
        
    #print matrice
    
    matricenewnop.sort()
    #matricenew.reverse()
    start_ord_operonsf=[]
    end_ord_operonsf=[]
    maxval_ord_operonsf=[]
    i=0
    for i in range(len(matricenewnop)):      
        start_ord_operonsf.append(matricenewnop[i][0])
        end_ord_operonsf.append(matricenewnop[i][1])
        maxval_ord_operonsf.append(matricenewnop[i][2])
    
    if len(maxval_ord_operonsf)>0:
        max_av=max(maxval_ord_operonsf)
        min_av=min(maxval_ord_operonsf)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Operons artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_operonsf.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_operonsf.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_operonsf)):
            medav=maxval_ord_operonsf[me]
            st=start_ord_operonsf[me]
            ed=end_ord_operonsf[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_operonsf.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_operonsf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_operonsf.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_operonsf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_operonsf.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_operonsf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_operonsf.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_operonsf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_operonsf.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_operonsf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_operonsf.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_operonsf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No operons"
    
    
    
    matricenewintf=[]
    pippo=0
    for pippo in range(len(start_interf)):
        array_internonintf=[]       
        array_internon0intf=start_interf[pippo]
        array_internon1intf=end_interf[pippo]
        array_internon2intf= max_interf[pippo]
        array_internonintf.append(array_internon0intf)
        array_internonintf.append(array_internon1intf)
        array_internonintf.append(array_internon2intf)
        matricenewintf.append(array_internonintf)
        
    #print matrice
    
    matricenewintf.sort()
    #matricenew.reverse()
    start_ord_interf=[]
    end_ord_interf=[]
    maxval_ord_interf=[]
    i=0
    for i in range(len(matricenewintf)):      
        start_ord_interf.append(matricenewintf[i][0])
        end_ord_interf.append(matricenewintf[i][1])
        maxval_ord_interf.append(matricenewintf[i][2])
    
    if len(maxval_ord_interf)>0:
        max_av=max(maxval_ord_interf)
        min_av=min(maxval_ord_interf)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Intergenic artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_interf.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_interf.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        for me in range(len(start_ord_interf)):
            medav=maxval_ord_interf[me]
            st=start_ord_interf[me]
            ed=end_ord_interf[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_interf.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_interf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_interf.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_interf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_interf.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_interf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_interf.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_interf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_interf.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_interf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_interf.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_interf.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No intergenic"
    
    
    
    
    
    
    
    
    
    
    
    matricenewr=[]
    pippo=0
    for pippo in range(len(start_utrr)):
        array_internonr=[]       
        array_internon0r=start_utrr[pippo]
        array_internon1r=end_utrr[pippo]
        array_internon2r= max_utrr[pippo]
        array_internonr.append(array_internon0r)
        array_internonr.append(array_internon1r)
        array_internonr.append(array_internon2r)
        matricenewr.append(array_internonr)
        
    #print matrice
    print "artemis files preparation R strand"
    matricenewr.sort()
    #matricenew.reverse()
    start_ord_utrr=[]
    end_ord_utrr=[]
    maxval_ord_utrr=[]
    i=0
    for i in range(len(matricenewr)):      
        start_ord_utrr.append(matricenewr[i][0])
        end_ord_utrr.append(matricenewr[i][1])
        maxval_ord_utrr.append(matricenewr[i][2])
    
    if len(maxval_ord_utrr)>0:
        max_av=max(maxval_ord_utrr)
        min_av=min(maxval_ord_utrr)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---UTR artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_utrr.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_utrr.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_utrr)):
            medav=maxval_ord_utrr[me]
            st=start_ord_utrr[me]
            ed=end_ord_utrr[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_utrr.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_utrr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_utrr.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_utrr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_utrr.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_utrr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_utrr.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_utrr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_utrr.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_utrr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_utrr.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_utrr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print"No UTR"
    
    
    
          
    matricer=[]
    pippo=0
    for pippo in range(len(start_utr_overr)):
        array_internor=[]        
        array_interno0r=start_utr_overr[pippo]
        array_interno1r=end_utr_overr[pippo]
        array_interno2r= max_utr_overr[pippo]
        array_internor.append(array_interno0r)
        array_internor.append(array_interno1r)
        array_internor.append(array_interno2r)
        matricer.append(array_internor)
        
    #print matrice
    
    matricer.sort()
    #matrice.reverse()
    start_ord_utr_overr=[]
    end_ord_utr_overr=[]
    maxval_ord_utr_overr=[]
    i=0
    for i in range(len(matricer)):      
        start_ord_utr_overr.append(matricer[i][0])
        end_ord_utr_overr.append(matricer[i][1])
        maxval_ord_utr_overr.append(matricer[i][2])
    
    if len(maxval_ord_utr_overr)>0:
        max_av=max(maxval_ord_utr_overr)
        min_av=min(maxval_ord_utr_overr)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Overlapping UTR artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_utr_overr.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_utr_overr.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_utr_overr)):
            medav=maxval_ord_utr_overr[me]
            st=start_ord_utr_overr[me]
            ed=end_ord_utr_overr[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_utr_overr.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_utr_overr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_utr_overr.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_utr_overr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_utr_overr.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_utr_overr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_utr_overr.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_utr_overr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_utr_overr.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_utr_overr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_utr_overr.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_utr_overr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        #for zz in range(len(start_utr_over)):
        #    file_utr_over.write("st %s  ---  end %s\n" %(start_utr_over[zz], end_utr_over[zz]))
    
    else:
        print "No overlapping UTR"
    
    matricenewnr=[]
    pippo=0
    for pippo in range(len(start_antir)):
        array_internonnr=[]       
        array_internonn0r=start_antir[pippo]
        array_internonn1r=end_antir[pippo]
        array_internonn2r= max_antir[pippo]
        array_internonnr.append(array_internonn0r)
        array_internonnr.append(array_internonn1r)
        array_internonnr.append(array_internonn2r)
        matricenewnr.append(array_internonnr)
        
    #print matrice
    
    matricenewnr.sort()
    #matricenew.reverse()
    start_ord_antir=[]
    end_ord_antir=[]
    maxval_ord_antir=[]
    i=0
    for i in range(len(matricenewnr)):      
        start_ord_antir.append(matricenewnr[i][0])
        end_ord_antir.append(matricenewnr[i][1])
        maxval_ord_antir.append(matricenewnr[i][2])
    
    if len(maxval_ord_antir)>0:
        max_av=max(maxval_ord_antir)
        min_av=min(maxval_ord_antir)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Antisense artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_antisenser.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_antisenser.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_antir)):
            medav=maxval_ord_antir[me]
            st=start_ord_antir[me]
            ed=end_ord_antir[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_antisenser.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_antisenser.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_antisenser.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_antisenser.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_antisenser.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_antisenser.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_antisenser.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_antisenser.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_antisenser.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_antisenser.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_antisenser.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_antisenser.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No antisense"
    
    
    
    
    
    matricenewnopr=[]
    pippo=0
    for pippo in range(len(start_operonsr)):
        array_internonnopr=[]       
        array_internonnop0r=start_operonsr[pippo]
        array_internonnop1r=end_operonsr[pippo]
        array_internonnop2r= max_operonsr[pippo]
        array_internonnopr.append(array_internonnop0r)
        array_internonnopr.append(array_internonnop1r)
        array_internonnopr.append(array_internonnop2r)
        matricenewnopr.append(array_internonnopr)
        
    #print matrice
    
    matricenewnopr.sort()
    #matricenew.reverse()
    start_ord_operonsr=[]
    end_ord_operonsr=[]
    maxval_ord_operonsr=[]
    i=0
    for i in range(len(matricenewnopr)):      
        start_ord_operonsr.append(matricenewnopr[i][0])
        end_ord_operonsr.append(matricenewnopr[i][1])
        maxval_ord_operonsr.append(matricenewnopr[i][2])
    
    if len(maxval_ord_operonsr)>0:
        max_av=max(maxval_ord_operonsr)
        min_av=min(maxval_ord_operonsr)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Operons artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_operonsr.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_operonsr.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        me=0
        for me in range(len(start_ord_operonsr)):
            medav=maxval_ord_operonsr[me]
            st=start_ord_operonsr[me]
            ed=end_ord_operonsr[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_operonsr.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_operonsr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_operonsr.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_operonsr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_operonsr.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_operonsr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_operonsr.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_operonsr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_operonsr.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_operonsr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_operonsr.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_operonsr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No operons"
    
    
    
    matricenewintr=[]
    pippo=0
    for pippo in range(len(start_interr)):
        array_internonintr=[]       
        array_internon0intr=start_interr[pippo]
        array_internon1intr=end_interr[pippo]
        array_internon2intr= max_interr[pippo]
        array_internonintr.append(array_internon0intr)
        array_internonintr.append(array_internon1intr)
        array_internonintr.append(array_internon2intr)
        matricenewintr.append(array_internonintr)
        
    #print matrice
    
    matricenewintr.sort()
    #matricenew.reverse()
    start_ord_interr=[]
    end_ord_interr=[]
    maxval_ord_interr=[]
    i=0
    for i in range(len(matricenewintr)):      
        start_ord_interr.append(matricenewintr[i][0])
        end_ord_interr.append(matricenewintr[i][1])
        maxval_ord_interr.append(matricenewintr[i][2])
    
    if len(maxval_ord_interr)>0:
        max_av=max(maxval_ord_interr)
        min_av=min(maxval_ord_interr)
        diff_av=max_av-min_av
        n_int=float(diff_av/6)
        print ("\n*---Intergenic artemis legend ---*")
        #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
        print ("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
        file_artemis_interr.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
        file_artemis_interr.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
        for me in range(len(start_ord_interr)):
            medav=maxval_ord_interr[me]
            st=start_ord_interr[me]
            ed=end_ord_interr[me]
            if medav >= min_av and medav <min_av+n_int:
                file_artemis_interr.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
                file_artemis_interr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int and medav <min_av+n_int+n_int:
                file_artemis_interr.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
                file_artemis_interr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
                file_artemis_interr.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
                file_artemis_interr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
                file_artemis_interr.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
                file_artemis_interr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
                file_artemis_interr.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
                file_artemis_interr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
            if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
                file_artemis_interr.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
                file_artemis_interr.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
    else:
        print "No intergenic"
    
    
    
    
    
    
    
    
    
    
    
    file_in1.close()
    file_in2.close()
    file_probes.close()
    
    file_utrf.close()
    file_utr_overf.close()
    file_antisensef.close()
    file_artemis_utrf.close()
    file_artemis_utr_overf.close()
    file_artemis_antisensef.close()
    
    file_utrr.close()
    file_utr_overr.close()
    file_antisenser.close()
    file_artemis_utrr.close()
    file_artemis_utr_overr.close()
    file_artemis_antisenser.close()
    

    print ("results stored into file '%s'. "
        % os.path.abspath(file_out_name1f) )
       


if __name__=='__main__':
    main()        
    