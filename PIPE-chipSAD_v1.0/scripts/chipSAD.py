import os
import sys
from optparse import OptionParser
from math import sqrt
from operator import itemgetter
from scipy import stats

def main():
    opts=parse_cmd_line()
    
    table(opts) 
    
    
    
def parse_cmd_line():
    #Parse the command line and return a parameters structure.

    usage = "usage: %prog [options] [structure]"
       
    parser = OptionParser(usage)
    parser.add_option("-o", "--outdir", action="store", dest="outdir", default=".", type="string", help="directory for output file. Default is current directory.")
    parser.add_option("-w", "--width", action="store", dest="width", default=40,  type="int", help="sliding window width. Default 40.")
    parser.add_option("-s", "--step", action="store", dest="step", default=10,  type="int", help="step of computation. Default 10.")
    parser.add_option("-t", "--threshold", action="store", dest="threshold", default=-300,  type="int", help="t-test threshold value. Default SAD compute it from ttest table.")
    parser.add_option("-p", "--percentile", action="store", dest="percentile", default=8,  type="int", help="percentile threshold. Default 8.")
    parser.add_option("-v", "--pvalue", action="store", dest="pvalue", default=0.05,  type="float", help="pvalue threshold. Default 0.05. If you set the t parameter different from defaul this options do not will be considered.")    
    parser.add_option("-g", "--gap", action="store", dest="gap", default=1,  type="int", help="the gap widht. Default 240.")
    parser.add_option("-b", "--pseudomedian", action="store", dest="pseudomedian", default=1,  type="int", help="boolean var: 0 for no presmoothing. Default 1.")
    parser.add_option("-m", "--minimum", action="store", dest="minimum", default=3,  type="int", help="minimun probe number for CPRs. Default 3.")
    parser.add_option("-j", "--strand", action="store", dest="strand", default="F",  type="string", help="strand. Default F.")



    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("An signal file, a probe file, the ttest table are required.")

    if not os.path.isdir(options.outdir):
        parser.error("Not a valid directory: '%s'." % options.outdir)

    options.entity = args[0]
    options.entity2= args[1]
    options.entity3=args[2]
    
    return options
    
 
def average(array):
    sum=0
    for h in range(len(array)):
        sum=sum+array[h]
    average=sum/len(array)
    return average

def st_dev(array, media):
    summ=0
    for o in range(len(array)):
        ter=(array[o]-media)**2
        summ=summ+ter       
    sdq=summ/len(array)
    sd=sqrt(sdq)
    return sd

def cov(array1, array2, media1, media2, Rcal):
    summ=0
    for oo in range(len(array1)):
        ter=((array1[oo]-media1)*(array2[oo]-media2))/((((array1[oo]-media1)**2)*((array2[oo]-media2)**2))+Rcal)
        summ=summ+ter   
    return summ

def correlazione (cov1, cov2):
    cor=(cov1*cov2)
    #/(sigma1*sigma2)
    return cor

def ttest_sw(media, sigma, mu, R):
    tsw=(media-mu)/(sigma+R)
    return tsw
    
def t_test(ave1, ave2, sigma1, sigma2, n, R):
    t=(ave2-ave1)/(sqrt((sigma1**2/n)+(sigma2**2/n))+R)
    return t
    
def t_test_nosigma(ave1, ave2, alpha0):
    t_ns=(ave2-ave1)/alpha0
    return t_ns

def den_def(sigma1, sigma2, n, R):
    den_noalpha=(sqrt((sigma1**2/n)+(sigma2**2/n)))+R
    return den_noalpha

def den_a(sigma1, sigma2, n):
    den_noalpha=(sqrt((sigma1**2/n)+(sigma2**2/n)))
    return den_noalpha


def alpha0_prova (sigma1, sigma2, n):
    a0=(sqrt((sigma1**2/n)+(sigma2**2/n)))
    return a0 
    
def num (ave1, ave2):
    num_ave=ave2-ave1
    return num_ave

def percentile(per, n):
    R=(per/100)*(n+1)
    return R

def per_cal(FR, eg, ep):
    Rcal=FR*(eg-ep)+ep
    return Rcal

def check(number):
    var=0
    if number%2==0:
        var=1
    return var

def pseudomedian (array_sw):
    pw_ave_array=[]
    for qq in range(len(array_sw)):
        elesw=array_sw[qq]
        for kk in range(len(array_sw)):
            if kk>=qq:
                pw_ave=float(elesw+array_sw[kk])/2
                pw_ave_array.append(pw_ave)
        
    pw_ave_array.sort()
    L=len(pw_ave_array)
    if check(L)==1:
        index=int(L/2)
        el1=pw_ave_array[index]
        el2=pw_ave_array[index+1]
        md=float(el1+el2)/2
    else:
        index=int(L/2)
        md=pw_ave_array[index]
    return md
 
 
def table(options):
    file_in_name = options.entity
    file_probes_name = options.entity2
    file_table_name = options.entity3

    file_out_base = os.path.splitext(os.path.basename(options.entity))[0]
    file_out_name = os.path.join(options.outdir,'ChipSAD_results_%s_strand%s_w%d_step%d_perc%d_th%d_pv%.3f_g%d_b%d_m%d.txt' % (file_out_base, options.strand, options.width, options.step, options.percentile, options.threshold, options.pvalue, options.gap, options.pseudomedian, options.minimum))
    file_ttest_name = os.path.join(options.outdir,'ChipSAD_ttest_%s_strand%s_w%d_step%d_perc%d_th%d_pv%.3f_g%d_b%d_m%d.dat' % (file_out_base, options.strand, options.width, options.step, options.percentile, options.threshold, options.pvalue, options.gap, options.pseudomedian, options.minimum))
    file_artemis_name = os.path.join(options.outdir,'ChipSAD_artemis_%s_strand%s_w%d_step%d_perc%d_th%d_pv%.3f_g%d_b%d_m%d.dat' % (file_out_base, options.strand, options.width, options.step, options.percentile, options.threshold, options.pvalue, options.gap, options.pseudomedian, options.minimum))                                                                                                                                                                                                                                                                  
    file_comparison_name = os.path.join(options.outdir,'ChipSAD_comparison_%s_strand%s_w%d_step%d_perc%d_th%d_pv%.3f_g%d_b%d_m%d.dat' % (file_out_base, options.strand, options.width, options.step, options.percentile, options.threshold, options.pvalue, options.gap, options.pseudomedian, options.minimum))                                                                                                                                                                                                                                                                  
    file_longprobe_name = os.path.join(options.outdir,'ChipSAD_probes_%s_strand%s_w%d_step%d_perc%d_th%d_pv%.3f_g%d_b%d_m%d.txt' % (file_out_base, options.strand, options.width, options.step, options.percentile, options.threshold, options.pvalue, options.gap, options.pseudomedian, options.minimum))                                                                                                                                                                                                                                                                  
    file_longprobe_name2 = os.path.join(options.outdir,'ChipSAD_excel_%s_strand%s_w%d_step%d_perc%d_th%d_pv%.3f_g%d_b%d_m%d.txt' % (file_out_base, options.strand, options.width, options.step, options.percentile, options.threshold, options.pvalue, options.gap, options.pseudomedian, options.minimum))                                                                                                                                                                                                                                                                  

    
    
    file_in = file(file_in_name, 'r')
    file_probes = file(file_probes_name, 'r')
    file_table= file(file_table_name, 'r')
    file_out = file(file_out_name, 'w')
    file_ttest = file(file_ttest_name, 'w')
    file_artemis = file(file_artemis_name, 'w')
    file_comparison=file(file_comparison_name, 'w')
    file_longprobe=file(file_longprobe_name, 'w')
    file_longprobe2=file(file_longprobe_name2, 'w')
    
    
    file_longprobe.write("strand\tCodingNonCoding\tORF\tpristart\tpristop\tReporterId\tM\tA\tTrascriptionalUnitNumber\tTU_start\tTU_end\tMeanValue\tTtest\tPvalue\n")
    
    
    mi_dat=[]
    Mvalues=[]
    array_average=[]
    array_stdev=[]
    array_average1=[]
    array_stdev1=[]
    array_len2=[]
    array_len1=[]
    
    count=0
    
    print ("\ndata storing")
    for line in file_in:
        mi_dat.append(float(line))
        Mvalues.append(float(line))
    
    while count<(options.width):
        ele=mi_dat[count]
        mi_dat.append(ele)
        count=count+1  
 
    
    value0_1_array=[]
    value0_05_array=[]
    value0_01_array=[]
    value0_001_array=[]
    gdl_table_array=[]
    
    for line in file_table:
        gdl_table=line.split()[0].split("\n")[0]
        value0_1=line.split()[1]
        value0_05=line.split()[2]
        value0_01=line.split()[3]
        value0_001=line.split()[4]
        gdl_table_array.append(int(gdl_table))
        value0_1_array.append(value0_1)
        value0_05_array.append(value0_05)
        value0_01_array.append(value0_01)
        value0_001_array.append(value0_001)


    
    start_probe_array=[]
    start_probe_sort_array=[]
    end_probe_array=[]
    reportedId_array=[]
    line_file=[]
    status_array=[]
    Mvalue_array=[]
    Avalue_array=[]
    kcont=0
    for line in file_probes:
        if options.strand == 'F':
            if (line.startswith('F')):
                sent=line.split()
                kconts= [sent.count(w) for w in sent]
                kcont=len(kconts)
                if kcont==8:
                    stp=int(line.split()[3])
                    endp=int(line.split()[4])
                    repId=line.split()[5]
                    status=line.split()[1]
                    Mvalue=float(line.split()[6])
                    Avalue=float(line.split()[7])
                    start_probe_array.append(stp)
                    end_probe_array.append(endp)
                    status_array.append(status)
                    Mvalue_array.append(Mvalue)
                    Avalue_array.append(Avalue)
                    reportedId_array.append(repId)
                    
                    l1=line.strip('\t')
                    l2=l1.strip('\n')
                    line_file.append(l2)
                
                if kcont==7:
                    stp=int(line.split()[2])
                    endp=int(line.split()[3])
                    repId=line.split()[4]
                    status=line.split()[1]
                    Mvalue=float(line.split()[5])
                    Avalue=float(line.split()[6])
                    start_probe_array.append(stp)
                    end_probe_array.append(endp)
                    reportedId_array.append(repId)
                    status_array.append(status)
                    Mvalue_array.append(Mvalue)
                    Avalue_array.append(Avalue)
                    l1=line.strip('\t')
                    l2=l1.strip('\n')
                    line_file.append(l2)
    
        if (options.strand == "R"):
            if (line.startswith('R')):
                sent=line.split()
                kconts= [sent.count(w) for w in sent]
                kcont=len(kconts)
                if kcont==8:
                    stp=int(line.split()[3])
                    endp=int(line.split()[4])
                    repId=line.split()[5]
                    status=line.split()[1]
                    Mvalue=float(line.split()[6])
                    Avalue=float(line.split()[7])
                    start_probe_array.append(stp)
                    end_probe_array.append(endp)
                    status_array.append(status)
                    Mvalue_array.append(Mvalue)
                    Avalue_array.append(Avalue)
                    reportedId_array.append(repId)
                    l1=line.strip('\t')
                    l2=l1.strip('\n')
                    line_file.append(l2)
                if kcont==7:
                    stp=int(line.split()[2])
                    endp=int(line.split()[3])
                    repId=line.split()[4]
                    status=line.split()[1]
                    Mvalue=float(line.split()[5])
                    Avalue=float(line.split()[6])
                    start_probe_array.append(stp)
                    end_probe_array.append(endp)
                    reportedId_array.append(repId)
                    status_array.append(status)
                    Mvalue_array.append(Mvalue)
                    Avalue_array.append(Avalue)
                    l1=line.strip('\t')
                    l2=l1.strip('\n')
                    line_file.append(l2)
   
   
    matrice=[]
    for pippo in range(len(start_probe_array)):
        array_interno=[]
        array_interno0= start_probe_array[pippo]
        array_interno1=end_probe_array[pippo]
        array_interno2=Mvalue_array[pippo]
        array_interno3=Avalue_array[pippo]
        array_interno4=line_file[pippo]
        array_interno.append(array_interno0)
        array_interno.append(array_interno1)
        array_interno.append(array_interno2)
        array_interno.append(array_interno3)
        array_interno.append(array_interno4)
        matrice.append(array_interno)
        
    
    
    matrice.sort()
    end_probe_sort_array=[]
    status_sort_array=[]
    start_probe_sort_array=[]
    Mvalue_sort_array=[]
    Avalue_sort_array=[]
    Mvalue_dasortare=[]
    line_file_sort_array=[]
    for i in range(len(matrice)):
        start_probe_sort_array.append(matrice[i][0])
        end_probe_sort_array.append(matrice[i][1])
        Mvalue_sort_array.append(matrice[i][2])
        Mvalue_dasortare.append(matrice[i][2])
        Avalue_sort_array.append(matrice[i][3])
        line_file_sort_array.append(matrice[i][4])
        
  
    
    
    
    avtot=average(Mvalue_sort_array)
    stdevtot=st_dev(Mvalue_sort_array, avtot)
    RcalM=stdevtot/10
    #print("RcalM %s\n" %(RcalM))
    
    
    
    
    tsw_th=0.3
    
    if options.pseudomedian == 1:
        start=0   
        th=options.width
        step=options.step
        #st=0
        
        print("CPRs selection and pseudomedian calculation")
        
        #while (st<= len(mi_dat)):  
               
        indicex=0
        prob=0    
        while(th<len(mi_dat)-1):
            #file_out.write( "entro while\n")
            #file_out.write("\n\nthreshold %d\n" %(th))
            #file_out.write("th %s, len %s\n" %(th, len(mi_dat)-options.width))
            array_sw=[] 
            array_sw1=[]
            array_swp=[] 
            array_sw1p=[]
            array_stp=[]
            array_endp=[]
            count_diffp=0
            array_stp2=[]
            array_endp2=[]
            count_diffp2=0
            tsw_array=[]
            tsw1_array=[]
            
            ks=0
            #if prob==1:
             #   z=0
            for z in range(indicex, len(start_probe_sort_array)-1):
                #file_out.write("pos %s\n" %(th))
                ddx=abs(start_probe_sort_array[z]-th)
                dsx=abs(start_probe_sort_array[z+1]-th)
                var=0
                #file_out.write("st_z %s, th %s, st_z+1 %s, ddx %s, dsx %s\n" %(start_probe_sort_array[z], th, start_probe_sort_array[z+1], ddx, dsx))
                
                if start_probe_sort_array[z]>=th and start_probe_sort_array[z+1]>=th:
                    #prob=start_probe_sort_array[z]-th
                    #th=th+step
                    #file_out.write("entro mag mag\n")
                    array_sw.append(-2)
                    array_sw1.append(-2)
                    array_sw.append(-2)
                    array_sw1.append(-2)
                    z=0
                    break
                    #var=1
                    
                
                if start_probe_sort_array[z]<=th and start_probe_sort_array[z+1]>=th and (ddx<=options.width or dsx<=options.width):
                #if start_probe_sort_array[z]<=th and start_probe_sort_array[z+1]>=th:
                    #file_out.write("z %s\n" %(z))
                    #file_out.write( "entro if\n")
                    prob=0
                    #file_out.write("st prima %s, st dopo %s\n" %(start_probe_sort_array[z], start_probe_sort_array[z+1]))
                    m=0
                    n=0
                    indicex=z
                    
                    #file_out.write("\nPRIMA: st %s end %s\n" %(start, th))
                    while count_diffp==0 :
                        ##probes a cavallo tra le due finestre
                        if start_probe_sort_array[z-m]<=th and end_probe_sort_array[z-m]>th:
                            #file_out.write("1 st %s, ed %s, th %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], th))
                            num_dx=th-start_probe_sort_array[z-m]
                            diff=end_probe_sort_array[z-m]-start_probe_sort_array[z-m]
                            I_dx=(num_dx*Mvalue_sort_array[z-m])/diff
                            #file_out.write("num_dx %s, diff %s, I_dx %s\n" %(num_dx, diff, I_dx))
                            array_swp.append(I_dx)
                            num_sx=end_probe_sort_array[z-m]-th
                            diff2=end_probe_sort_array[z-m]-start_probe_sort_array[z-m]
                            I_sx=(num_sx*Mvalue_sort_array[z-m])/diff2
                            #file_out.write("num_sx %s, diff %s, I_sx %s\n" %(num_sx, diff2, I_sx))
                            array_sw1p.append(I_sx)
                            
                            if len(array_swp)>options.minimum:
                                last=array_swp[len(array_swp)-1]
                                array_swp.pop(len(array_swp)-1)
                                media=average(array_swp)
                                st_sw=st_dev(array_swp, media)                               
                                if last<=media+(1*st_sw) and last>=media-(1*st_sw):
                                    array_sw.append(I_dx)
                                    array_swp.append(I_dx)
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(th)
                                    #file_out.write("1 cavallo selez st %s, ed %s, last %s\n" %(start_probe_sort_array[z-m], th, last))
                                    m=m+1
                                else:
                                    count_diffp=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_swp)-1))
                                
                            else:
                                array_sw.append(I_dx)
                                array_stp.append(start_probe_sort_array[z-m])
                                array_endp.append(th)
                                #file_out.write("prima obb cavallo 1: st %s, ed %s, I_dx %s\n" %(start_probe_sort_array[z-m], th, I_dx))
                                m=m+1
                                
                            if len(array_sw1p)>options.minimum:
                                last1=array_sw1p[len(array_sw1p)-1]
                                array_sw1p.pop(len(array_sw1p)-1)
                                media1=average(array_sw1p)                               
                                st_sw1=st_dev(array_sw1p, media1)
                                #tsw1=ttest_sw(media1, st_sw1,  I_sx, RcalM)
                                #tsw1_array.append(tsw1)
                                ##file_out.write("2 cavallo candidata st %s, ed %s, tsw %s, len %s\n" %(th, end_probe_sort_array[z-m],  tsw, len(array_sw1p)))
                                #gdl=len(array_sw1p)-1
                                #file_out.write("2 candidata cavallo st %s, ed %s, last %s, len %s\n" %(th, end_probe_sort_array[z-m],  last1, len(array_sw1p)))
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                
                                #if tsw1_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw1<=tsw_th and tsw1>=tsw_th_neg:
                                #if tsw1>=tsw_th:
                                    array_stp2.append(th)
                                    array_endp2.append(end_probe_sort_array[z-m])
                                    array_sw1.append(I_sx)
                                    #file_out.write("2 cavallo selez st %s, ed %s, tsw %s\n" %(th, end_probe_sort_array[z-m],  tsw))
                                    m=m+1
                                else:
                                    count_diffp2=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_sw1p)-1))
                                """
                                if last1<=media1+(1*st_sw1) and last1>=media1-(1*st_sw1):
                                    array_sw1.append(I_sx)
                                    array_sw1p.append(I_sx)
                                    array_endp2.append(end_probe_sort_array[z-m])
                                    array_sw1.append(I_sx)
                                    #file_out.write("2 cavallo selez st %s, ed %s, last %s\n" %(th, end_probe_sort_array[z-m], last1))
                                    m=m+1
                                else:
                                    count_diffp2=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_sw1p)-1))
                                    
                            else:
                                array_stp2.append(th)
                                array_endp2.append(end_probe_sort_array[z-m])
                                array_sw1.append(I_sx)
                                #file_out.write("cavallo obb 2: st %s, ed %s, I_sx %s\n" %(th, end_probe_sort_array[z-m],  I_sx))
                                m=m+1
                            
                        else:
                            ##probes non a cavallo
                            #file_out.write("1 st %s, ed %s, th %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], th))
                            array_swp.append(Mvalue_sort_array[z-m])
                            if len(array_swp)>options.minimum:
                                last=array_swp[len(array_swp)-1]
                                array_swp.pop(len(array_swp)-1)
                                media=average(array_swp)
                                st_sw=st_dev(array_swp, media)
                                #tsw=ttest_sw(media, st_sw,  Mvalue_sort_array[z-m], RcalM)
                                #tsw_array.append(tsw)
                                #file_out.write("1 candidata st %s, ed %s, last %s, len %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], last, len(array_swp)))
                                #gdl=len(array_swp)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                
                                #if tsw_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw<=tsw_th and tsw>=tsw_th_neg:
                                    #if tsw>=tsw_th:
                                    #file_out.write("1 selz st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], tsw))
                                    array_sw.append(Mvalue_sort_array[z-m])
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(end_probe_sort_array[z-m])
                                    m=m+1
                                else:
                                    count_diffp=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_swp)-1))
                                """
                                if last<=media+(1*st_sw) and last>=media-(1*st_sw):
                                    array_sw.append(Mvalue_sort_array[z-m])
                                    array_swp.append(Mvalue_sort_array[z-m])
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(end_probe_sort_array[z-m])
                                    #file_out.write("1 selz st %s, ed %s, last %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], last))
                                    m=m+1
                                else:
                                    count_diffp=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_swp)-1))
                                
                                
                            else:
                                #file_out.write("selz obb 1: st %s, ed %s, M %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], Mvalue_sort_array[z-m]))
                                array_sw.append(Mvalue_sort_array[z-m])
                                array_stp.append(start_probe_sort_array[z-m])
                                array_endp.append(end_probe_sort_array[z-m])
                                m=m+1
                            
                    #file_out.write("SECONDA: st %s end %s\n" %(start+options.width, th+options.width))     
                    ##seconda finestra
                    while count_diffp2==0 and z+n<len(start_probe_sort_array):
                        #file_out.write("2 %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                        array_sw1p.append(Mvalue_sort_array[z+n])
                        if len(array_sw1p)>options.minimum:
                            last1=array_sw1p[len(array_sw1p)-1]
                            array_sw1p.pop(len(array_sw1p)-1)
                            media1=average(array_sw1p)
                            st_sw1=st_dev(array_sw1p, media1)
                            #tsw1=ttest_sw(media1, st_sw1,  Mvalue_sort_array[z+n], RcalM)
                            #tsw1_array.append(tsw1)
                            #file_out.write("2 cadidata st %s, ed %s, last %s, len %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], last1, len(array_sw1p)))
                            #gdl=len(array_sw1p)-1
                            """
                            for val in range(len(gdl_table_array)):
                                #print "entro2"
                                if gdl==gdl_table_array[val]:
                                    #print "entro3"
                                    index_gdl=gdl_table_array.index(gdl_table_array[val])
                                    tsw_th=float(value0_05_array[index_gdl])
                                    #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                if gdl>=121:
                                    index_gdl=34
                                    tsw_th=float(value0_05_array[index_gdl])
                            
                            #if tsw1_array[0]>0:
                            tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                            if tsw1<=tsw_th and tsw1>=tsw_th_neg:
                            #if tsw1>=tsw_th:
                                #file_out.write("selz 2 st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], tsw))
                                array_stp2.append(start_probe_sort_array[z+n])
                                array_endp2.append(end_probe_sort_array[z+n])
                                array_sw1.append(Mvalue_sort_array[z+n])
                                n=n+1
                            else:
                                count_diffp2=1
                                #file_out.write("nella finestra ci sono %s probes\n" %(len(array_sw1p)-1))
                            """
                            if last1<=media1+(1*st_sw1) and last1>=media1-(1*st_sw1):
                                #file_out.write("selz 2 st %s, ed %s, last %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], last1))
                                array_stp2.append(start_probe_sort_array[z+n])
                                array_endp2.append(end_probe_sort_array[z+n])
                                array_sw1.append(Mvalue_sort_array[z+n])
                                array_sw1p.append(Mvalue_sort_array[z+n])
                                n=n+1
                            else:
                                count_diffp2=1
                                #file_out.write("nella finestra ci sono %s probes\n" %(len(array_sw1p)-1))
                        else:
                            #file_out.write("obb 2: st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                            array_stp2.append(start_probe_sort_array[z+n])
                            array_endp2.append(end_probe_sort_array[z+n])
                            array_sw1.append(Mvalue_sort_array[z+n])
                            n=n+1
                    break
                ##file_out.write ("len %s\n" %(len(start_probe_sort_array)-2))
                
                
                #break
                
               
                if start_probe_sort_array[z]<=th and start_probe_sort_array[z+1]>=th and (ddx>options.width and dsx>options.width):
                    #file_out.write("\nELSE PRIMA: st %s end %s\n" %(start, th))
                    #file_out.write("ELSE SECONDA: st %s end %s\n" %(start+options.width, th+options.width))  
                    array_sw.append(-2)
                    array_sw1.append(-2)
                    array_sw.append(-2)
                    array_sw1.append(-2)
                    var=1
                    z=z+1
                    break
                
                
                
                    
                
                    
                if z == len(start_probe_sort_array)-2 and var==0:
                #if z == len(start_probe_sort_array)-2:
                    #file_out.write("\nultima\n")
                    #print "entro if 2"
                    if th >=start_probe_sort_array[z] and th >=start_probe_sort_array[z+1] and start_probe_sort_array[0]<=options.width: 
                        #print "entro else"
                        n=0
                        m=1
                        
                        #file_out.write("PRIMA: st %s end %s\n" %(start, th))
                        while count_diffp==0 :
                            #file_out.write(" 1 st %s, ed %s, th %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], th))
                            array_swp.append(Mvalue_sort_array[z-m])
                            if len(array_swp)>options.minimum:
                                last=array_swp[len(array_swp)-1]
                                array_swp.pop(len(array_swp)-1)
                                media=average(array_swp)
                                st_sw=st_dev(array_swp, media)
                                #tsw=ttest_sw(media, st_sw,  Mvalue_sort_array[z-m], RcalM)
                                #tsw_array.append(tsw)
                                #file_out.write("1 candidata st %s, ed %s, last %s, len %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], last, len(array_swp)))
                                #gdl=len(array_swp)-1
                                """
                                
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                
                                #if tsw_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw<=tsw_th and tsw>=tsw_th_neg:
                                    #if tsw>=tsw_th:
                                    #file_out.write("1 sel st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], tsw))
                                    array_sw.append(Mvalue_sort_array[z-m])
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(end_probe_sort_array[z-m])
                                    m=m+1
                                else:
                                    count_diffp=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_swp)-1))
                                """
                                if last<=media+(1*st_sw) and last>=media-(1*st_sw):
                                    #file_out.write("1 sel st %s, ed %s, last %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], last))
                                    array_sw.append(Mvalue_sort_array[z-m])
                                    array_swp.append(Mvalue_sort_array[z-m])
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(end_probe_sort_array[z-m])
                                    m=m+1
                                else:
                                    count_diffp=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_swp)-1))
                            else:
                                #file_out.write("obb 1: st %s, ed %s, last %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], last))
                                array_sw.append(Mvalue_sort_array[z-m])
                                array_stp.append(start_probe_sort_array[z-m])
                                array_endp.append(end_probe_sort_array[z-m])
                                m=m+1
                            
                            """
                            ##file_out.write( "entrowh\n")
                            array_stp.append(start_probe_sort_array[len(start_probe_sort_array)-m])
                            #print "pippo"
                            array_endp.append(end_probe_sort_array[len(end_probe_sort_array)-m])
                            array_sw.append(Mvalue_sort_array[len(Mvalue_sort_array)-m])
                            diff=end_probe_sort_array[len(start_probe_sort_array)-m]-start_probe_sort_array[len(end_probe_sort_array)-m]
                            count_diffp=count_diffp+diff
                            m=m+1
                            #print "pippo"
                            ##file_out.write("stp %s, edp %s, diff %s, count %s, m %s\n" %(start_probe_sort_array[len(start_probe_sort_array)-m], end_probe_sort_array[len(end_probe_sort_array)-m], diff, count_diffp,  m))
                            """
                        #file_out.write("SECONDA: st %s end %s\n" %(start+options.width, th+options.width))
                        z=0
                        while count_diffp2==0 and z+n<len(start_probe_sort_array):
                            #file_out.write("2 st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                            array_sw1p.append(Mvalue_sort_array[z+n])
                            if len(array_sw1p)>options.minimum:
                                last1=array_sw1p[len(array_sw1p)-1]
                                array_sw1p.pop(len(array_sw1p)-1)
                                media1=average(array_sw1p)
                                st_sw1=st_dev(array_sw1p, media1)
                                #tsw1=ttest_sw(media1, st_sw1,  Mvalue_sort_array[z+n], RcalM)
                                #tsw1_array.append(tsw1)
                                #file_out.write("2 candidata st %s, ed %s, last %s, len %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], last1, len(array_sw1p)))
                                #gdl=len(array_sw1p)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                
                                #if tsw1_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw1<=tsw_th and tsw1>=tsw_th_neg:
                                    #if tsw1>=tsw_th:
                                    #file_out.write("2 sel st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], tsw))
                                    array_stp2.append(start_probe_sort_array[z+n])
                                    array_endp2.append(end_probe_sort_array[z+n])
                                    array_sw1.append(Mvalue_sort_array[z+n])
                                    n=n+1
                                else:
                                    count_diffp2=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_sw1p)-1))
                                """
                                if last1<=media1+(1*st_sw1) and last1>=media1-(1*st_sw1):
                                    #file_out.write("2 sel st %s, ed %s, last %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], last1))
                                    array_stp2.append(start_probe_sort_array[z+n])
                                    array_endp2.append(end_probe_sort_array[z+n])
                                    array_sw1.append(Mvalue_sort_array[z+n])
                                    array_sw1p.append(Mvalue_sort_array[z+n])
                                    n=n+1
                                else:
                                    count_diffp2=1
                                    #file_out.write("nella finestra ci sono %s probes\n" %(len(array_sw1p)-1))
                                
                            else:
                                #file_out.write("2 obb: st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                                array_stp2.append(start_probe_sort_array[z+n])
                                array_endp2.append(end_probe_sort_array[z+n])
                                array_sw1.append(Mvalue_sort_array[z+n])
                                n=n+1
                        break
                    
                    else:
                        array_sw.append(-2)
                        array_sw1.append(-2)
                        array_sw.append(-2)
                        array_sw1.append(-2)
                            
                        """
                        array_stp2.append(start_probe_sort_array[z+n])
                        array_endp2.append(end_probe_sort_array[z+n])
                        array_sw1.append(Mvalue_sort_array[z+n])
                        diff2=end_probe_sort_array[z+n]-start_probe_sort_array[z+n]
                        count_diffp2=count_diffp2+diff2
                        n=n+1
                        """
                
                
            
            
            """           
            #for j in range(start-enlarge,th):
                ##file_out.write("PRIMA W: j %d, start %d, th %d\n" %(j, start-enlarge, th))
            for s_e in range(len(array_stp)):
                for con in range(array_stp[s_e], array_endp[s_e]):
                    elem=mi_dat[con]
                    array_sw.append(elem)
                    #print "entro"
                    ##file_out.write("ele %f\n" %(elem))
                ##file_out.write("st prima window %s, ed prima window %s\n" %(array_stp[s_e], array_endp[s_e]) )
                
            for s_e2 in range(len(array_stp2)):
                for con2 in range(array_stp2[s_e2], array_endp2[s_e2]):
                    elem2=mi_dat[con2]
                    array_sw1.append(elem2)
                    ##file_out.write("ele %f\n" %(elem2))
                ##file_out.write("st seconda window %s, ed seconda window %s\n" %(array_stp2[s_e2], array_endp2[s_e2]) )
               
                #elem1=mi_dat[j+options.width+enlarge2]
                #array_sw1.append(elem1)   
            """
            
            ave=pseudomedian(array_sw)
            ave1=pseudomedian(array_sw1)
            stdev=st_dev(array_sw, ave)
            stdev1=st_dev(array_sw1, ave1)
            array_average.append(ave)
            array_stdev.append(stdev)
            array_average1.append(ave1)
            array_stdev1.append(stdev1)
            array_len1.append(len(array_sw))
            array_len2.append(len(array_sw1))
            #file_out.write("pseudomedian %f, stdev %f\n" %(ave, stdev))
            #file_out.write("pseudomedian1 %f, stdev1 %f\n" %(ave1, stdev1))
            start=start+step
            th=th+step
            #st=st+step
                    #file_out.write("start %d, th %d\n" %(start, th))
                #file_out.write("count %d, elem %f\n" %(count, elem) )
        

    
    
    else:
        ###questa parte va ancora modificata### 
        start=0   
        th=options.width
        step=options.step
        #st=0
        enlarge=0
        enlarge2=0
        indice=0
        indice2=0
        wind_st=0
        wind_ed=800
        indicex=0
        print("average and standard deviation calculation")
        
        #while (st<= len(mi_dat)):  
        ks=0       
        while(th<len(mi_dat)-1):
            #file_out.write( "entro while\n")
            #file_out.write("\n\nthreshold %d\n" %(th))
            #file_out.write("th %s, len %s\n" %(th, len(mi_dat)-options.width))
            array_sw=[] 
            array_sw1=[]
            array_stp=[]
            array_endp=[]
            count_diffp=0
            array_stp2=[]
            array_endp2=[]
            count_diffp2=0
            ks=0
            tsw_array=[]
            tsw1_array=[]
            array_swp=[]
            array_sw1p=[]
            ks=0
            
            for z in range(indicex, len(start_probe_sort_array)-1):
                #file_out.write("pos %s\n" %(th))
                if start_probe_sort_array[z]<=th and start_probe_sort_array[z+1]>=th:
                    #file_out.write("z %s\n" %(z))
                    #file_out.write( "entro if\n")
                    #file_out.write("st prima %s, st dopo %s\n" %(start_probe_sort_array[z], start_probe_sort_array[z+1]))
                    m=0
                    n=0
                    indicex=z
                    
                    #file_out.write("\nPRIMA: st %s end %s\n" %(start, th))
                    while count_diffp==0 :
                        if start_probe_sort_array[z-m]<=th and end_probe_sort_array[z-m]>th:
                            #file_out.write("st %s, ed %s, th %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], th))
                            num_dx=th-start_probe_sort_array[z-m]
                            diff=end_probe_sort_array[z-m]-start_probe_sort_array[z-m]
                            I_dx=(num_dx*Mvalue_sort_array[z-m])/diff
                            #file_out.write("num_dx %s, diff %s, I_dx %s\n" %(num_dx, diff, I_dx))
                            array_swp.append(I_dx)
                            num_sx=end_probe_sort_array[z-m]-th
                            diff2=end_probe_sort_array[z-m]-start_probe_sort_array[z-m]
                            I_sx=(num_sx*Mvalue_sort_array[z-m])/diff2
                            #file_out.write("num_sx %s, diff %s, I_sx %s\n" %(num_sx, diff2, I_sx))
                            array_sw1p.append(I_sx)
                            
                            if len(array_swp)>1:
                                media=average(array_swp)
                                st_sw=st_dev(array_swp, media)                               
                                tsw=ttest_sw(media, st_dev(array_swp, media),  I_dx, RcalM)
                                tsw_array.append(tsw)
                                #file_out.write("media %s, stdev %s, tsw %s\n" %(media, st_sw, tsw))
                                gdl=len(array_swp)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                """
                                #if tsw_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw<=tsw_th and tsw>=tsw_th_neg:
                                    #if tsw>=tsw_th:
                                    array_sw.append(I_dx)
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(th)
                                    #file_out.write("1 st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], th, tsw))
                                    m=m+1
                                else:
                                    count_diffp=1
                                
                            else:
                                array_sw.append(I_dx)
                                array_stp.append(start_probe_sort_array[z-m])
                                array_endp.append(th)
                                #file_out.write("1 len < 0: st %s, ed %s, I_dx %s\n" %(start_probe_sort_array[z-m], th, I_dx))
                                m=m+1
                                
                            if len(array_sw1p)>1:
                                media1=average(array_sw1p)                               
                                st_sw1=st_dev(array_sw1p, media1)
                                tsw1=ttest_sw(media1, st_sw1,  I_sx, RcalM)
                                tsw1_array.append(tsw1)
                                #file_out.write("media1 %s, stdev1 %s, tsw1 %s\n" %(media1, st_sw1, tsw1))
                                gdl=len(array_sw1p)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                """
                                #if tsw1_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw1<=tsw_th and tsw1>=tsw_th_neg:
                                    #if tsw1>=tsw_th:
                                    array_stp2.append(th)
                                    array_endp2.append(end_probe_sort_array[z-m])
                                    array_sw1.append(I_sx)
                                    #file_out.write("2 st %s, ed %s, tsw %s\n" %(th, end_probe_sort_array[z-m],  tsw,))
                                    m=m+1
                                else:
                                    count_diffp2=1
                                
                                    
                            else:
                                array_stp2.append(th)
                                array_endp2.append(end_probe_sort_array[z-m])
                                array_sw1.append(I_sx)
                                #file_out.write("2 len < 0: st %s, ed %s, I_sx %s\n" %(th, end_probe_sort_array[z-m],  I_sx))
                                m=m+1
                            
                        else:
                            #file_out.write("1 else st %s, ed %s, th %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], th))
                            array_swp.append(Mvalue_sort_array[z-m])
                            if len(array_swp)>1:
                                media=average(array_swp)
                                st_sw=st_dev(array_swp, media)
                                tsw=ttest_sw(media, st_sw,  Mvalue_sort_array[z-m], RcalM)
                                tsw_array.append(tsw)
                                #file_out.write("media %s, stdev %s, tsw %s\n" %(media, st_sw, tsw))
                                gdl=len(array_swp)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                """
                                #if tsw_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw<=tsw_th and tsw>=tsw_th_neg:
                                    #if tsw>=tsw_th:
                                    #file_out.write("1 else st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], tsw))
                                    array_sw.append(Mvalue_sort_array[z-m])
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(end_probe_sort_array[z-m])
                                    m=m+1
                                else:
                                    count_diffp=1
                                
                            else:
                                #file_out.write("1 else len < 0: st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], tsw))
                                array_sw.append(Mvalue_sort_array[z-m])
                                array_stp.append(start_probe_sort_array[z-m])
                                array_endp.append(end_probe_sort_array[z-m])
                                m=m+1
                            
                    #file_out.write("SECONDA: st %s end %s\n" %(start+options.width, th+options.width))     
                    while count_diffp2==0 and z+n<len(start_probe_sort_array):
                        #file_out.write("2 st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                        array_sw1p.append(Mvalue_sort_array[z+n])
                        if len(array_sw1p)>1:
                            media1=average(array_sw1p)
                            st_sw1=st_dev(array_sw1p, media1)
                            tsw1=ttest_sw(media1, st_sw1,  Mvalue_sort_array[z+n], RcalM)
                            tsw1_array.append(tsw1)
                            #file_out.write("media1 %s, stdev1 %s, tsw1 %s\n" %(media1, st_sw1, tsw1))
                            gdl=len(array_sw1p)-1
                            """
                            for val in range(len(gdl_table_array)):
                                #print "entro2"
                                if gdl==gdl_table_array[val]:
                                    #print "entro3"
                                    index_gdl=gdl_table_array.index(gdl_table_array[val])
                                    tsw_th=float(value0_05_array[index_gdl])
                                    #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                if gdl>=121:
                                    index_gdl=34
                                    tsw_th=float(value0_05_array[index_gdl])
                            """
                            #if tsw1_array[0]>0:
                            tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                            if tsw1<=tsw_th and tsw1>=tsw_th_neg:
                            #if tsw1>=tsw_th:
                                #file_out.write("2 st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], tsw))
                                array_stp2.append(start_probe_sort_array[z+n])
                                array_endp2.append(end_probe_sort_array[z+n])
                                array_sw1.append(Mvalue_sort_array[z+n])
                                n=n+1
                            else:
                                count_diffp2=1
                            
                        else:
                            #file_out.write("2 len < 0 : st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                            array_stp2.append(start_probe_sort_array[z+n])
                            array_endp2.append(end_probe_sort_array[z+n])
                            array_sw1.append(Mvalue_sort_array[z+n])
                            n=n+1
                    break
                #file_out.write ("len %s\n" %(len(start_probe_sort_array)-2))
                
                if z == len(start_probe_sort_array)-2:
                    #file_out.write("\nultima\n")
                    #print "entro if 2"
                    if th >=start_probe_sort_array[z] and th >=start_probe_sort_array[z+1]: 
                        #print "entro else"
                        n=0
                        m=1
                        
                        #file_out.write("PRIMA: st %s end %s\n" %(start, th))
                        while count_diffp==0 :
                            #file_out.write("1 else st %s, ed %s, th %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], th))
                            array_swp.append(Mvalue_sort_array[z-m])
                            if len(array_swp)>1:
                                media=average(array_swp)
                                st_sw=st_dev(array_swp, media)
                                tsw=ttest_sw(media, st_sw,  Mvalue_sort_array[z-m], RcalM)
                                tsw_array.append(tsw)
                                #file_out.write("media %s, stdev %s, tsw %s\n" %(media, st_sw, tsw))
                                gdl=len(array_swp)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                """
                                #if tsw_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw<=tsw_th and tsw>=tsw_th_neg:
                                    #if tsw>=tsw_th:
                                    #file_out.write("1 else st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], tsw))
                                    array_sw.append(Mvalue_sort_array[z-m])
                                    array_stp.append(start_probe_sort_array[z-m])
                                    array_endp.append(end_probe_sort_array[z-m])
                                    m=m+1
                                else:
                                    count_diffp=1
                                
                            else:
                                #file_out.write("1 else len < 0: st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z-m], end_probe_sort_array[z-m], tsw))
                                array_sw.append(Mvalue_sort_array[z-m])
                                array_stp.append(start_probe_sort_array[z-m])
                                array_endp.append(end_probe_sort_array[z-m])
                                m=m+1
                            
                            """
                            #file_out.write( "entrowh\n")
                            array_stp.append(start_probe_sort_array[len(start_probe_sort_array)-m])
                            #print "pippo"
                            array_endp.append(end_probe_sort_array[len(end_probe_sort_array)-m])
                            array_sw.append(Mvalue_sort_array[len(Mvalue_sort_array)-m])
                            diff=end_probe_sort_array[len(start_probe_sort_array)-m]-start_probe_sort_array[len(end_probe_sort_array)-m]
                            count_diffp=count_diffp+diff
                            m=m+1
                            #print "pippo"
                            #file_out.write("stp %s, edp %s, diff %s, count %s, m %s\n" %(start_probe_sort_array[len(start_probe_sort_array)-m], end_probe_sort_array[len(end_probe_sort_array)-m], diff, count_diffp,  m))
                            """
                        #file_out.write("SECONDA: st %s end %s\n" %(start+options.width, th+options.width))
                        z=0
                        while count_diffp2==0 and z+n<len(start_probe_sort_array):
                            #file_out.write("2 st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                            array_sw1p.append(Mvalue_sort_array[z+n])
                            if len(array_sw1p)>1:
                                media1=average(array_sw1p)
                                st_sw1=st_dev(array_sw1p, media1)
                                tsw1=ttest_sw(media1, st_sw1,  Mvalue_sort_array[z+n], RcalM)
                                tsw1_array.append(tsw1)
                                #file_out.write("media1 %s, stdev1 %s, tsw1 %s\n" %(media1, st_sw1, tsw1))
                                gdl=len(array_sw1p)-1
                                """
                                for val in range(len(gdl_table_array)):
                                    #print "entro2"
                                    if gdl==gdl_table_array[val]:
                                        #print "entro3"
                                        index_gdl=gdl_table_array.index(gdl_table_array[val])
                                        tsw_th=float(value0_05_array[index_gdl])
                                        #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                                    if gdl>=121:
                                        index_gdl=34
                                        tsw_th=float(value0_05_array[index_gdl])
                                """
                                #if tsw1_array[0]>0:
                                tsw_th_neg=float(tsw_th)-(float(tsw_th)+float(tsw_th))
                                if tsw1<=tsw_th and tsw1>=tsw_th_neg:
                                    #if tsw1>=tsw_th:
                                    #file_out.write("2 st %s, ed %s, tsw %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], tsw))
                                    array_stp2.append(start_probe_sort_array[z+n])
                                    array_endp2.append(end_probe_sort_array[z+n])
                                    array_sw1.append(Mvalue_sort_array[z+n])
                                    n=n+1
                                else:
                                    count_diffp2=1
                                
                            else:
                                #file_out.write("2 len < 0 : st %s, ed %s, th %s\n" %(start_probe_sort_array[z+n], end_probe_sort_array[z+n], th))
                                array_stp2.append(start_probe_sort_array[z+n])
                                array_endp2.append(end_probe_sort_array[z+n])
                                array_sw1.append(Mvalue_sort_array[z+n])
                                n=n+1
                        break
                            
                            
                        """
                        array_stp2.append(start_probe_sort_array[z+n])
                        array_endp2.append(end_probe_sort_array[z+n])
                        array_sw1.append(Mvalue_sort_array[z+n])
                        diff2=end_probe_sort_array[z+n]-start_probe_sort_array[z+n]
                        count_diffp2=count_diffp2+diff2
                        n=n+1
                        """
                
                
            
            
            """            
            #for j in range(start-enlarge,th):
                #file_out.write("PRIMA W: j %d, start %d, th %d\n" %(j, start-enlarge, th))
            for s_e in range(len(array_stp)):
                for con in range(array_stp[s_e], array_endp[s_e]):
                    elem=mi_dat[con]
                    array_sw.append(elem)
                    #print "entro"
                    #file_out.write("ele %f\n" %(elem))
                #file_out.write("st prima window %s, ed prima window %s\n" %(array_stp[s_e], array_endp[s_e]) )
                
            for s_e2 in range(len(array_stp2)):
                for con2 in range(array_stp2[s_e2], array_endp2[s_e2]):
                    elem2=mi_dat[con2]
                    array_sw1.append(elem2)
                    #file_out.write("ele %f\n" %(elem2))
                #file_out.write("st seconda window %s, ed seconda window %s\n" %(array_stp2[s_e2], array_endp2[s_e2]) )
               
                #elem1=mi_dat[j+options.width+enlarge2]
                #array_sw1.append(elem1)
            
            for hh in range(len(array_sw)):
                file_out.write("1 st %s, end %s\n" %(array_stp[hh], array_endp[hh]))
            for hhh in range(len(array_sw1)):
                file_out.write("2 st %s, end %s\n" %(array_stp2[hhh], array_endp2[hhh]))
            """
            ave=average(array_sw)
            ave1=average(array_sw1)
            stdev=st_dev(array_sw, ave)
            stdev1=st_dev(array_sw1, ave1)
            array_average.append(ave)
            array_stdev.append(stdev)
            array_average1.append(ave1)
            array_stdev1.append(stdev1)
            array_len1.append(len(array_sw))
            array_len2.append(len(array_sw1))
            #file_log.write("average %f, stdev %f\n" %(ave, stdev))
            #file_log.write("average1 %f, stdev1 %f\n" %(ave1, stdev1))
            start=start+step
            th=th+step
            #st=st+step
            #file_out.write("start %d, th %d\n" %(start, th))
            #file_out.write("ks %s\n" %(ks))
            #wind_st=wind_st+(indice2-ks)
            #if (wind_ed+indice2 <= len(start_probe_sort_array)-1):
             #   wind_ed=wind_ed+ks
            #else:
                #print "finito"
             #   wind_ed = len(start_probe_sort_array)
                #wind_st=wind_st+indice-ks
            """
            if start in range(options.width, len(start_probe_array), options.width ):
                wind_st=wind_st+1
                if (wind_ed+options.width <= len(start_probe_sort_array)-1):
                    wind_ed=wind_ed+options.width
                    #file_out.write("count %d, elem %f\n" %(count, elem) )
                else:
                    wind_ed=len(start_probe_sort_array)-1
            """
    #print len(array_average)
    #print len(array_average1)
    #print len(array_stdev)
    #print len(array_stdev1)
    
    count_pos=0
    count_neg=0 
    count_explosion=0
    new_array_stdev=[]
    array_a0=[]
    ttest_array=[]
    num_array=[]
    den_array=[]
    
    
    for dd in range(len(array_stdev)):
        if (dd != len(array_stdev)):
            if (array_stdev[dd]!=0 and array_stdev1[dd]!=0):
                a0_p=den_a(array_stdev[dd], array_stdev1[dd], options.width)
                #a0_p=alpha0_prova(array_stdev[dd], array_stdev1[dd], options.width)
                if a0_p>0:
                    array_a0.append(a0_p) 
                    #file_log.write("elem %f\n" %(a0_p))
   
    array_a0.sort()
    #print len(array_a0)
    #print ("percentile calculation")
    R=percentile(float(options.percentile), float(len(array_a0)))
    RI=int( R)
    if RI<R:
        Rsx=RI+1
        Rdx=RI
        FR=R-RI
    else: 
        Rsx=RI
        Rdx=RI-1
        FR=RI-R

    #file_log.write("R %f, Rsx %d, Rdx %d, FR %f\n " %(R, Rsx, Rdx, FR))
     
    #for v in range(len(array_a0)):
    ele_g=array_a0[Rsx]
    ele_p=array_a0[Rdx]  
    Rcal=per_cal(FR, ele_g, ele_p)
    #file_log.write("ele_g %.20f, ele_p %.20f, Rcal %.20f\n " %(ele_g, ele_p, Rcal))
        
    

    min_a0=(min(array_a0))
    #file_log.write("min a0 %.50f\n" %(min_a0))
    
    for g in range(len(array_average)):
        if (g != len(array_average)):
            n_ave=num(array_average[g], array_average1[g])
            for p in range (0, step):
                num_array.append(n_ave)
                #file_num.write("%.10f\n" %(n_ave))
    #file_num.write("lungh norm %d\n" %(len(num_array)))
    #file_log.write("\nmin_num %f, max_num %f\n" %(min(num_array), max(num_array)))  
    #print len(num_array)
    
    print("ttest computation")
    len_array=[]
    for z in range(len(array_average)):
        if (z != len(array_average)):
            #file_out.write("av %s, av1 %s\n" %(array_average[z], array_average1[z]))
            if array_average[z] == -2 and array_average1[z]== -2:
                ttest=0
            else:
                ttest=t_test(array_average[z], array_average1[z], array_stdev[z], array_stdev1[z], options.width, Rcal)
            den=den_def(array_stdev[z], array_stdev1[z], options.width, Rcal)
            lm=(array_len1[z]+array_len2[z])-2
            for k in range (0, step):
                #file_out.write("%.10f\n" %(ttest))  
                ttest_array.append(ttest)
                #file_den.write("%.10f\n" %(den))
                den_array.append(den)
                len_array.append(lm)
           
    
    #file_log.write("\nmin_ttest %f, max_ttest %f\n" %(min(ttest_array), max(ttest_array)))   
    #file_log.write("\nmin_den %f, max_den %f\n" %(min(den_array), max(den_array)))  
    
    ##file_out.write("\n\ncount_pos %d, count_neg %d, count_explosion %d\n" %(count_pos, count_neg, count_explosion))  
    #file_log.write("\n\ncount_pos %d, count_neg %d\n" %(count_pos, count_neg))          
    #file_ttest.write("\ncount ttest %d\n" %(len(ttest_array)))  
    #file_den.write("\ncount den %d\n" %(len(den_array)))
    
    #print len(den_array)
    
    ttest_array_new=[]
    num_array_new=[]
    den_array_new=[]
    len_array_new=[]
    st=len(ttest_array)-options.width
    for r in range(st,len(ttest_array)):
        ttest_array_new.append(ttest_array[r])
        num_array_new.append(num_array[r])
        den_array_new.append(den_array[r])
        len_array_new.append(len_array[r])
        
    for rr in range(0,st):
        ttest_array_new.append(ttest_array[rr])
        num_array_new.append(num_array[rr])
        den_array_new.append(den_array[rr])
        len_array_new.append(len_array[rr])
    
    #file_log.write("len ttest %d, len ttest_new %d\n" %(len(ttest_array), len(ttest_array_new)))
    
    for rrr in range(len(ttest_array_new)):
        file_ttest.write("%f\n" %(ttest_array_new[rrr]))
        #file_num.write("%f\n" %(num_array_new[rrr]))
        #file_den.write("%f\n" %(den_array_new[rrr]))  
    
    #print len(ttest_array_new)

    
    
        #print gdl_table
    
    
    
    
    print("choise of signal areas ")
    k_pos=0 
    k_start=0
    start_array=[]
    end_array=[]
    M_media_array=[]
    M_min_array=[]
    M_max_array=[]
    M_sigma_array=[]
    M_cov_array=[]
    matrice_cor=[]
    ttest_th=0
    validita=0
    a=0
    countator=0
    indice_st=0
    for u in range(len(ttest_array_new)):
        if (u != len(ttest_array_new)-1):
            
            gdl=len_array_new[u]
            #print ("ttest %s, gdl array %s\n" %(ttest_array_new[u], gdl))
            
            if options.threshold == -300:
                
                if options.pvalue == 0.1:
                    for val in range(len(gdl_table_array)):
                        if gdl==gdl_table_array[val]:
                            index_gdl=gdl_table_array.index(gdl_table_array[val])
                            ttest_th=float(value0_1_array[index_gdl])
                        if gdl>=121:
                            index_gdl=34
                            ttest_th=float(value0_1_array[index_gdl])
                        
                if options.pvalue == 0.05:
                    #print "entro1"
                    for val in range(len(gdl_table_array)):
                        #print "entro2"
                        if gdl==gdl_table_array[val]:
                            #print "entro3"
                            index_gdl=gdl_table_array.index(gdl_table_array[val])
                            ttest_th=float(value0_05_array[index_gdl])
                            #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                        if gdl>=121:
                            index_gdl=34
                            ttest_th=float(value0_05_array[index_gdl])
                            #print ("gdl %s, ttest_th %s\n" %(gdl, ttest_th))
                            
                        
                if options.pvalue == 0.01:
                    for val in range(len(gdl_table_array)):
                        if gdl==gdl_table_array[val]:
                            index_gdl=gdl_table_array.index(gdl_table_array[val])
                            ttest_th=float(value0_01_array[index_gdl])
                        if gdl>=121:
                            index_gdl=34
                            ttest_th=float(value0_01_array[index_gdl])
                        
                if options.pvalue == 0.001:
                    for val in range(len(gdl_table_array)):
                        if gdl==gdl_table_array[val]:
                            index_gdl=gdl_table_array.index(gdl_table_array[val])
                            ttest_th=float(value0_001_array[index_gdl])
                        if gdl>=121:
                            index_gdl=34
                            ttest_th=float(value0_001_array[index_gdl])
                        
            else:
                ttest_th=float(options.threshold)
     
            
            ttest_th_neg=float(ttest_th)-(float(ttest_th)+float(ttest_th))
            
            if ttest_array_new[u]>=0:
            
                if ttest_array_new[u]>=ttest_th:
                    k_start=k_start+1
            
                
                #if ttest_array_new[u]>ttest_th_neg and ttest_array_new[u]<ttest_th:
                 #   k_start=k_start+1
                if ttest_array_new[u]<=ttest_th and k_start>0:
                    #k_pos=k_pos+1
                    array_temp_M=[]
                    end_elem=ttest_array_new[u]
                    end_index=u+2  
                    start_elem=ttest_array_new[u-k_start-1]
                    start_index=u-k_start+1
                    #for ind in (start_index, end_index):
                    #file_out.write("\nstart %d  end  %d, length %d\n" %(start_index, end_index, k_start))
                    for a in range(indice_st, len(start_probe_sort_array)):
                        #file_out.write("st %s, end %s,    start_p %d  end_p  %d, M %s\n" %(start_index, end_index, start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))

                        if (start_probe_sort_array[a]>=start_index and end_probe_sort_array[a]<=end_index) or ( (start_probe_sort_array[a]<=start_index and end_probe_sort_array[a]>= start_index) or (start_probe_sort_array[a]<=end_index and end_probe_sort_array[a]>=end_index)): 
                            #file_out.write("st %s, end %s,    start_p %d  end_p  %d, M %s\n" %(start_index, end_index, start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))

                            if ( (start_probe_sort_array[a]<start_index and (end_probe_sort_array[a]-start_index)>=(((end_probe_sort_array[a]-start_probe_sort_array[a])*30)/100)) or (end_probe_sort_array[a]>end_index and ((end_index-start_probe_sort_array[a])>=(((end_probe_sort_array[a]-start_probe_sort_array[a])*30)/100)))):               
                                array_temp_M.append(Mvalue_sort_array[a])
                                #file_out.write("start_p %d  end_p  %d, M %s\n" %(start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                            if (start_probe_sort_array[a]>=start_index) and end_probe_sort_array[a]<=end_index:
                                array_temp_M.append(Mvalue_sort_array[a])
                                #file_out.write("start_p %d  end_p  %d, M %s\n" %(start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                            #file_out.write("start_p %d  end_p  %d, M %s\n" %(start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                        #indice_st=a
                        #break
                    if len(array_temp_M)!=0:
                        media_M=average(array_temp_M)
                        #file_out.write("media %s\n" %(media_M))
                        #cov_M=cov(array_temp_M, media_M)
                        max_M=max(array_temp_M)
                        min_M=min(array_temp_M)
                        sigma_M=st_dev(array_temp_M, media_M)
                        M_media_array.append(media_M)
                        M_sigma_array.append(sigma_M)
                            #M_cov_array.append(cov_M)
                        M_max_array.append(max_M)
                        M_min_array.append(min_M)
                        end_array.append(end_index)
                        start_array.append(start_index)
                        for el in range(len(array_temp_M)):
                            array_int=[]
                            array_interno=array_temp_M[el]
                            array_int.append(array_interno)
                            matrice_cor.append(array_int)
                        #file_out.write("END: start %d  end  %d, media %s, sigma %s, n val %s\n" %(start_index, end_index, media_M, sigma_M, len(array_temp_M)))
                        k_start=0
            
            
            if ttest_array_new[u]<0:
            
                if ttest_array_new[u]<=ttest_th_neg:
                    k_start=k_start+1
            
                #ttest_th_neg=float(ttest_th)-(float(ttest_th)+float(ttest_th))
                #if ttest_array_new[u]>ttest_th_neg and ttest_array_new[u]<ttest_th:
                 #   k_start=k_start+1
                if ttest_array_new[u]>ttest_th_neg and k_start>0:
                    #k_pos=k_pos+1
                    array_temp_M=[]
                    end_elem=ttest_array_new[u]
                    end_index=u+2  
                    start_elem=ttest_array_new[u-k_start-1]
                    start_index=u-k_start+1
                    #for ind in (start_index, end_index):
                    #file_out.write("\nstart %d  end  %d, length %d\n" %(start_index, end_index, k_start))
                    for a in range(indice_st, len(start_probe_sort_array)):
                        #file_out.write("st %s, end %s,    start_p %d  end_p  %d, M %s\n" %(start_index, end_index, start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                            
                        if (start_probe_sort_array[a]>=start_index and end_probe_sort_array[a]<=end_index) or ( (start_probe_sort_array[a]<=start_index and end_probe_sort_array[a]>= start_index) or (start_probe_sort_array[a]<=end_index and end_probe_sort_array[a]>=end_index)): 
                            #file_out.write("st %s, end %s,    start_p %d  end_p  %d, M %s\n" %(start_index, end_index, start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                            #array_temp_M.append(Mvalue_sort_array[a])
                            #"""
                            if ( (start_probe_sort_array[a]<start_index and (end_probe_sort_array[a]-start_index)>=(((end_probe_sort_array[a]-start_probe_sort_array[a])*30)/100)) or (end_probe_sort_array[a]>end_index and ((end_index-start_probe_sort_array[a])>=(((end_probe_sort_array[a]-start_probe_sort_array[a])*30)/100)))):               
                                array_temp_M.append(Mvalue_sort_array[a])
                                #file_out.write("start_p %d  end_p  %d, M %s\n" %(start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                            if (start_probe_sort_array[a]>=start_index) and end_probe_sort_array[a]<=end_index:
                                array_temp_M.append(Mvalue_sort_array[a])
                                #file_out.write("start_p %d  end_p  %d, M %s\n" %(start_probe_sort_array[a], end_probe_sort_array[a], Mvalue_sort_array[a]))
                            #"""
                        #indice_st=a
                        #break
                    if len(array_temp_M)!=0:
                        media_M=average(array_temp_M)
                        #file_out.write("media %s\n" %(media_M))
                        #cov_M=cov(array_temp_M, media_M)
                        max_M=max(array_temp_M)
                        min_M=min(array_temp_M)
                        sigma_M=st_dev(array_temp_M, media_M)
                        M_media_array.append(media_M)
                        M_sigma_array.append(sigma_M)
                            #M_cov_array.append(cov_M)
                        M_max_array.append(max_M)
                        M_min_array.append(min_M)
                        end_array.append(end_index)
                        start_array.append(start_index)
                        for el in range(len(array_temp_M)):
                            array_int=[]
                            array_interno=array_temp_M[el]
                            array_int.append(array_interno)
                            matrice_cor.append(array_int)
                        #file_out.write("END: start %d  end  %d, media %s, sigma %s, n val %s\n" %(start_index, end_index, media_M, sigma_M, len(array_temp_M)))
                        k_start=0
           
    #print len(start_array)
    #print len(end_array)
    start_array.append(start_array[len(start_array)-1])
    end_array.append(end_array[len(end_array)-1])
    M_media_array.append(M_media_array[len(M_media_array)-1])
    M_sigma_array.append(M_sigma_array[len(M_sigma_array)-1])
    M_max_array.append(M_max_array[len(M_max_array)-1])
    M_min_array.append(M_min_array[len(M_min_array)-1])
    
    #print len(start_array)
    
    print("artemis file calculation")
    k_diff=0
    array_max=[]
    ttest_av=[]
    ttest_new_array=[]
    start_new_array=[]
    end_new_array=[]
    ttest_dict={}
    v=0
    array1=[]
    array2=[]
    array_temp_m_tot=[]
    array_temp_m=[]
    temp_m_v=[]
    temp_m_v1=[]
    for v in range(len(start_array)):
        
        if(v != len(start_array)-1):
            #print v
            #array1.append(matrice_cor[v][0])
            #array2.append(matrice_cor[v+1][0])
            diff=start_array[v+1]-end_array[v]
            #file_out.write("\nstart+1 %d end %d diff %d\n" %( start_array[v+1], end_array[v], diff))
            #m_ave=abs((M_media_array[v+1]+M_media_array[v])/2)
            
            if M_media_array[v]>=0:
                ele_array_v=M_max_array[v]
            if M_media_array[v]<0:
                ele_array_v=M_min_array[v]
            if M_media_array[v+1]>=0:
                ele_array_v1=M_max_array[v+1]
            if M_media_array[v+1]<0:
                ele_array_v1=M_min_array[v+1]
            m_ave_tot=abs((ele_array_v1+ele_array_v)/2)
            m_diff=abs(ele_array_v1-ele_array_v)
            temp_m_v.append(ele_array_v)
            temp_m_v1.append(ele_array_v1)
           
            
           
            if(diff<=options.gap) and v != len(start_array)-2:
                #print 'gap'
                k_diff=k_diff+1
                #file_out.write("v1>v, opt, kdiff %d\n" %(k_diff))
                                        
            #file_out.write("kdiff %d\n" %(k_diff))
            else:
                #print 'entro'
                #if (M_media_array[v+1]-(M_sigma_array[v+1]*2))<=(M_media_array[v]+(M_sigma_array[v]*2)):
                #if m_diff< m_ave_tot:
                m_ave_tot_pu=abs((temp_m_v1[len(temp_m_v1)-1]+temp_m_v[0])/2)
                m_diff_pu=abs(temp_m_v1[len(temp_m_v1)-1]-temp_m_v[0])
                #file_out.write("m_diff %s, m_ave %s, m_diff_pu %s, m_ave_pu %s\n" %(m_diff, m_ave_tot, m_diff_pu, m_ave_tot_pu))
                media_temp=average(temp_m_v)
                if media_temp>=0:
                    picco=max(temp_m_v)
                else:
                    picco=min(temp_m_v)
                #if m_diff< m_ave_tot and m_diff_pu < m_ave_tot_pu and m_diff<2 and m_diff_pu<2:
                if m_diff< m_ave_tot and m_diff_pu < m_ave_tot_pu and abs(picco/2)<abs(temp_m_v1[len(temp_m_v1)-1]) and v != len(start_array)-2:
                    #file_out.write("M_media_v %s,  M_media_v1 %s\n" %(M_media_array[v], M_media_array[v+1] ))
                    #file_out.write("m_diff %s, m_ave %s\n" %(m_diff, m_ave_tot))
                    k_diff=k_diff+1
                    #file_out.write("kdiff %d\n" %(k_diff))
                else:
                    #file_out.write("else M_media_v %s,  M_media_v1 %s\n" %(M_media_array[v], M_media_array[v+1] ))
                    #file_out.write("else m_diff %s, m_ave %s\n" %(m_diff, m_ave_tot))
                    #file_out.write(" else, kdiff %d\n" %(k_diff))
                    #print "entro"
                    #k_diff=0
                    st_new=v-k_diff
                    start_index_new=start_array[st_new]
                    end_index_new=end_array[v]
                    #print start_index_new
                    start_new_array.append(start_index_new)
                    end_new_array.append(end_index_new)
                    ##file_out.write("start %d end %d\n" %(start_index_new, end_index_new))
                    #file_out.write("st_new %s, start %d end %d" %(st_new, start_index_new, end_index_new))
                    
                    for w in range(start_index_new, end_index_new):
                        #file_out.write("ttest elem %f\n" %(ttest_array_new[w]))
                        #ttest_av.append(float(ttest_array_new[w]))
                        if Mvalues[w]!=0.0:
                            ttest_av.append(float(Mvalues[w]))
                    ttest_av_new=pseudomedian(ttest_av)
                    ttest_new_array.append(ttest_av_new)
                    
                    k_diff=0
                    ttest_av=[]
                    array_temp_m=[]
                    array_temp_m_tot=[]
                    temp_m_v=[]
                    temp_m_v1=[]
    
    
    
    stm=0
    stp=0
    var=0
    Sval=0
    #print len(Mvalues)
    for Mval in range (stm, len(Mvalues)):
        #file_comparison.write("Mval %s\n" %(Mval))
        
        #if Sval==len(start_new_array) and Mval<len(Mvalues):
         #   file_comparison.write("0\n")
        for Sval in range(stp, len(start_new_array)):
            #file_comparison.write("start %s\n" %(start_new_array[Sval]))
            if Mval==start_new_array[Sval]:
                #file_comparison.write("end %s\n" %(end_new_array[Sval]))
                for pippo in range(start_new_array[Sval], end_new_array[Sval]):
                    file_comparison.write("%s\n" %(ttest_new_array[Sval]))
                    stm=stm+1
                var=end_new_array[Sval]
                stp=Sval+1
            
            if Mval>=var and Mval!=start_new_array[Sval]:
                file_comparison.write("0\n")
            break
    
    
    #print len(end_new_array)
    
    for topolino in range(end_new_array[len(end_new_array)-1], len(Mvalues)):
        file_comparison.write("0\n")
    #for pippo in range(len(array_max)):
     #   file_comparison.write("%s\n" %(array_max[pippo]))
            
    #print len(start_new_array)
    #print len(end_new_array)
    #print len(ttest_new_array)
    
         
    max_av=max(ttest_new_array)
    min_av=min(ttest_new_array)
    diff_av=max_av-min_av
    n_int=float(diff_av/6)
    print ("\n*--- artemis legend ---*")
    #print ("max %f, min %f, diff %f, n_int %f\n " %(max_av, min_av, diff_av, n_int))
    print("*--- blue<%.3f, light blue<%.3f, green<%.3f, yellow<%.3f, orange<%.3f, red<%.3f, black->end ---*\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
    #file_artemis.write("# BASE VAL1 VAL2 VAL3 VAL4 VAL5 VAL6 VAL7\n")
    file_artemis.write("# BASE VAL1<%.3f VAL2<%.3f VAL3<%.3f VAL4<%.3f VAL5<%.3f VAL6<%.3f VAL7end\n" %(min_av+n_int, min_av+n_int+n_int, min_av+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int, min_av+n_int+n_int+n_int+n_int+n_int, max_av))
    file_artemis.write("# colour 0:0:255 0:255:255 0:255:0 255:255:0 244:164:96 255:0:0 100:100:100 \n")  
    for me in range(len(start_new_array)):
        medav=ttest_new_array[me]
        st=start_new_array[me]
        ed=end_new_array[me]
        if medav >= min_av and medav <min_av+n_int:
            file_artemis.write("%d %f 0 0 0 0 0 0\n" %(st, medav))
            file_artemis.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        if medav >= min_av+n_int and medav <min_av+n_int+n_int:
            file_artemis.write("%d 0 %f 0 0 0 0 0\n" %(st, medav))
            file_artemis.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        if medav >= min_av+n_int+n_int and medav <min_av+n_int+n_int+n_int:
            file_artemis.write("%d 0 0 %f 0 0 0 0\n" %(st, medav))
            file_artemis.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        if medav >= min_av+n_int+n_int+n_int and medav <min_av+n_int+n_int+n_int+n_int:
            file_artemis.write("%d 0 0 0 %f 0 0 0\n" %(st, medav))
            file_artemis.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        if medav >= min_av+n_int+n_int+n_int+n_int and medav < min_av+n_int+n_int+n_int+n_int+n_int:
            file_artemis.write("%d 0 0 0 0 %f 0 0\n" %(st, medav))
            file_artemis.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        if medav >= min_av+n_int+n_int+n_int+n_int+n_int and medav <= max_av:
            file_artemis.write("%d 0 0 0 0 0 %f 0\n" %(st, medav))
            file_artemis.write("%d 0 0 0 0 0 0 %f\n" %(ed, medav))
        
    
    
    matricenew=[]
    pippo=0
    for pippo in range(len(ttest_new_array)):
        array_internon=[]
        array_internon0= ttest_new_array[pippo]
        array_internon1=start_new_array[pippo]
        array_internon2=end_new_array[pippo]
        array_internon.append(array_internon0)
        array_internon.append(array_internon1)
        array_internon.append(array_internon2)
        matricenew.append(array_internon)
        
    #print matrice
    
    matricenew.sort()
    matricenew.reverse()
    start_ord_array=[]
    end_ord_array=[]
    ttest_ord_array=[]
    i=0
    for i in range(len(matricenew)):
        ttest_ord_array.append(matricenew[i][0])
        start_ord_array.append(matricenew[i][1])
        end_ord_array.append(matricenew[i][2])
        
        
    
    
    
    
    
    print("results files compilation\n")
    probes_dict={}
     
    for m in range(len(reportedId_array)):
        rId_ele=reportedId_array[m]
        lf_ele=line_file[m]
        probes_dict[rId_ele]=[lf_ele]
    
    
    tu_name=0
    start_probe_tu=[]
    end_probe_tu=[]
    nome_probe_tu=[]
    tuname=[]
    value_tu=[]
    ttest_tu=[]
    pvalue_tu=[]
    start_array_tu=[]
    end_array_tu=[]
    for q in range(len(start_ord_array)):
        countN=0
        countC=0
        Mvalue_array_C=[]
        Mvalue_array_N=[]
        Avalue_array_C=[]
        Avalue_array_N=[]
        Mval_array_tot=[]
        Aval_array_tot=[]
        start=start_ord_array[q]
        end=end_ord_array[q]
        ele=ttest_ord_array[q]
        array_tm=[]
        array_str=[]
        tu_name_var=''
        tu_name=tu_name+1
        file_out.write("\n\n*** start %d end %d max value, average, t-test, p-value %f " %(start, end, ele))
        for u in range(len(start_probe_array)):
            
            #int_stp=int(start_probe_array[u])
            #int_endp=int(end_probe_array[u])
            #file_out.write("start prob %d end prob %d\n" %(start_probe_array[u],  end_probe_array[u]))
            #file_out.write("start prob %d end prob %d\n" %(int_stp,  int_endp))
            if (start_probe_array[u]>=start_ord_array[q] and end_probe_array[u]<=end_ord_array[q]) or ( (start_probe_array[u]<=start_ord_array[q] and end_probe_array[u]>= start_ord_array[q]) or (start_probe_array[u]<=end_ord_array[q] and end_probe_array[u]>=end_ord_array[q])):
                #file_out.write("start prob %d end prob %d\n" %(start_probe_array[u],  end_probe_array[u]))
                #rId=reportedId_array[u]
                #string=str(probes_dict[rId])
                if ( (start_probe_array[u]<start_ord_array[q] and (end_probe_array[u]-start_ord_array[q])>=(((end_probe_array[u]-start_probe_array[u])*30)/100)) or (end_probe_array[u]>end_ord_array[q] and ((end_ord_array[q]-start_probe_array[u])>=(((end_probe_array[u]-start_probe_array[u])*30)/100)))):               
                    array_tm.append(Mvalue_array[u])
            
            #if (start_probe_array[u]>=start_ord_array[q] and end_probe_array[u]<=end_ord_array[q]) or ( (start_probe_array[u]<=start_ord_array[q] and end_probe_array[u]>= start_ord_array[q]) or (start_probe_array[u]<=end_ord_array[q] and end_probe_array[u]>=end_ord_array[q])):

                    string=str(line_file[u])
                    array_str.append(string)
                    #file_out.write("%s\n" %(string))
                    if status_array[u]== 'C':
                        countC=countC+1
                        Mvalue_array_C.append(Mvalue_array[u])
                        Avalue_array_C.append(Avalue_array[u])
                        #file_out.write("elem array %f" %(Mvalue_array[u]))
                    if status_array[u]== 'N':
                        countN=countN+1
                        Mvalue_array_N.append(Mvalue_array[u])
                        Avalue_array_N.append(Avalue_array[u])
                                    
                if (start_probe_array[u]>=start_ord_array[q]) and end_probe_array[u]<=end_ord_array[q]:
                    
                
                
                    array_tm.append(Mvalue_array[u])
                
            #if     (start_probe_array[u]>=start_ord_array[q] and end_probe_array[u]<=end_ord_array[q]) or ( (start_probe_array[u]<=start_ord_array[q] and end_probe_array[u]>= start_ord_array[q]) or (start_probe_array[u]<=end_ord_array[q] and end_probe_array[u]>=end_ord_array[q])):

                    string=str(line_file[u])
                    array_str.append(string)
                    #file_out.write("%s\n" %(string))
                    if status_array[u]== 'C':
                        countC=countC+1
                        Mvalue_array_C.append(Mvalue_array[u])
                        Avalue_array_C.append(Avalue_array[u])
                        #file_out.write("elem array %f" %(Mvalue_array[u]))
                    if status_array[u]== 'N':
                        countN=countN+1
                        Mvalue_array_N.append(Mvalue_array[u])
                        Avalue_array_N.append(Avalue_array[u])
        Mval_array_tot.extend(Mvalue_array_C)
        Mval_array_tot.extend(Mvalue_array_N)
        Aval_array_tot.extend(Avalue_array_C)
        Aval_array_tot.extend(Avalue_array_N)
        av_m=average(array_tm)
        file_out.write("%s " %(av_m))
        file_out.write("%6.3f %6.4f\n" %  stats.ttest_1samp(array_tm, 0.0))
        for uffa in array_str:
            file_out.write("%s\n" %(uffa))
            start_tu=uffa.split("\t")[3]
            end_tu=uffa.split("\t")[4]
            nome_tu=uffa.split("\t")[5]
            tu_name_var='SAS_'+str(tu_name)+'_'+options.strand
            file_longprobe.write("%s\t%s\t%s\t%s\t%s\t" %(uffa, tu_name_var, start, end, av_m))
            file_longprobe.write("%6.3f\t%e\n" %  stats.ttest_1samp(array_tm, 0.0))
            start_probe_tu.append(int(start_tu))
            end_probe_tu.append(int(end_tu))
            nome_probe_tu.append(nome_tu)
            tuname.append(tu_name_var)
            value_tu.append(av_m)
            ttest, pvalue=stats.ttest_1samp(array_tm, 0.0)
            ttest_tu.append(ttest)
            pvalue_tu.append(pvalue)
            start_array_tu.append(start)
            end_array_tu.append(end)
        if countC>0 and countN>0:
            file_out.write("\nttest (MC+MN) = %6.3f  pvalue (MC+MN) = %6.4f\n" %  stats.ttest_1samp(Mval_array_tot, 0.0))
            file_out.write("ttest (AC+AN) = %6.3f  pvalue (AC+AN) = %6.4f\n" %  stats.ttest_1samp(Aval_array_tot, 0.0))
        if countC>0:
            aveC=average(Mvalue_array_C)
            aveCA=average(Avalue_array_C)
            file_out.write("\ncount C = %d, Maverage C = %.3f, Aaverage C = %.3f\n" %(countC, aveC, aveCA))
            #if countN>0:
            file_out.write("ttest MC = %6.3f  pvalue MC = %6.4f\n" %  stats.ttest_1samp(Mvalue_array_C, 0.0))
            file_out.write("ttest AC = %6.3f  pvalue AC = %6.4f\n" %  stats.ttest_1samp(Avalue_array_C, 0.0))
        else:
            file_out.write("\ncount C = 0\n")
        if countN>0:
            aveN=average(Mvalue_array_N)
            aveNA=average(Avalue_array_N)
            file_out.write("\ncount N = %d, Maverage N = %.3f, Aaverage N = %.3f\n" %(countN, aveN, aveNA))
            #if countC>0:
            file_out.write("ttest MN = %6.3f  pvalue MN = %6.4f\n" %  stats.ttest_1samp(Mvalue_array_N, 0.0))
            file_out.write("ttest AN = %6.3f  pvalue AN = %6.4f\n" %  stats.ttest_1samp(Avalue_array_N, 0.0))
        else:
            file_out.write("\ncount N = 0\n")
    
    
    
    file_longprobe2.write("strand\tstart_probe\tend_probe\tSAS_name\tSAS_start\tSAS_end\tSAS_pseudomedian_value\tSAS_ttest_value\tSAS_pvalue\n")
    i=0
    j=0
    var=0
    var2=0
    for i in range(len(start_probe_array)):
        var=0
        var2=0
        for j in range(len(start_probe_tu)):
            if start_probe_array[i]==start_probe_tu[j] and var2==0:
                file_longprobe2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(options.strand, start_probe_array[i], end_probe_array[i], reportedId_array[i], tuname[j], start_array_tu[j], end_array_tu[j], value_tu[j], ttest_tu[j], pvalue_tu[j]))
                var=1
                var2=1
        if var==0:
            file_longprobe2.write("%s\t%s\t%s\t%s\tNO_SAS\n" %(options.strand, start_probe_array[i], end_probe_array[i], reportedId_array[i]))
    
    
    
 

    file_in.close()
    #file_out.close()
    #file_num.close()
    file_ttest.close()
    file_artemis.close()
    #file_log.close()

    print ("results stored into file '%s'. "
        % os.path.abspath(file_out_name) )
       


if __name__=='__main__':
    main()        