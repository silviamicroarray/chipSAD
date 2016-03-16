import os
import sys
from optparse import OptionParser
from math import sqrt
from operator import itemgetter
#from scipy import stats
import fnmatch
import igraph
import igraph.clustering



def main():
    opts=parse_cmd_line()
    
    table(opts) 
    

def average(array):
    sum=0
    for h in range(len(array)):
        sum=sum+array[h]
    average=sum/len(array)
    return average

def make_list(size):
    """create a list of size number of zeros"""
    mylist = []
    for i in range(size):
        mylist.append(0)
    return mylist
 
def make_matrix(rows, cols):
    """
    create a 2D matrix as a list of rows number of lists
    where the lists are cols in size
    resulting matrix contains zeros
    """
    matrix = []
    for i in range(rows):
        matrix.append(make_list(cols))
    return matrix

def average_mod(st_end_array, mval_array):
    Mval_sum=sum(mval_array)
    N=len(st_end_array)
    peso=[]
    i=0
    for i in range(len(mval_array)):
        peso.append(float(mval_array[i])/float(Mval_sum))
    #print peso
    j=0
    den=0
    for j in range(len(peso)):
        den=den+(float(st_end_array[j])*peso[j])
        #print den
    media=float(den)/float(sum(peso))
    return media 
    
def parse_cmd_line():
    #Parse the command line and return a parameters structure.

    usage = "usage: %prog [options] [structure]"
       
    parser = OptionParser(usage)
    parser.add_option("-o", "--outdir", action="store", dest="outdir", default=".", type="string", help="directory for output file. Default is current directory.")
    parser.add_option("-u", "--thup", action="store", dest="thup", default=0.00, type="float", help="up threshold for DE tramscripts. Default is 0.00.")
    parser.add_option("-d", "--thdown", action="store", dest="thdown", default=-0.00, type="float", help="down threshold for DE tramscripts. Default is 0.00.")
    parser.add_option("-s", "--strand", action="store", dest="strand", default="F", type="string", help="set the strand. default F.")


    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("files to be analysed should be in the folder and not as argument of command")
    if not os.path.isdir(options.outdir):
        parser.error("Not a valid directory: '%s'." % options.outdir)

    #options.entity1 = args[0]
    
   
    
    return options
    
    
    
def table(options):
    #file_in_name1 = os.path.join(options.outdir, options.entity1)
    
  
    

    #file_out_base = os.path.splitext(os.path.basename(options.entity1))[0]
    #file_out_base2 = os.path.splitext(os.path.basename(options.entity2))[0]
    
    file_out_name=os.path.join(options.outdir,'%s_results_alignedTU.txt' %(options.strand) )
    file_report_name=os.path.join(options.outdir,'%s_report_alignedTU.txt' %(options.strand) )
    
    #file_in1 = file(file_in_name1, 'r')
    
    
    file_out = file(file_out_name, 'w')
    file_report = file(file_report_name, 'w') 
    
    file_report.write("Exp_name\tNum_up_regolted\tNum_down_regolated\t\n")
    count_exp=0
    strand_array=[]
    nameinter_array=[]
    orfup_array=[]
    orfdown_array=[]
    schema_array=[]
    start_array=[]
    end_array=[]
    length_array=[]
    Mvalue_array=[]
    ttest_array=[]
    pvalue_array=[]
    

    for filename in os.listdir(os.getcwd()):
    	file_in = file(filename, 'r')
        count_up=0
        count_down=0
        count_exp=count_exp+1
        #print filename
        line_array=[]
        for line in file_in:
            ln=(line.split("\n"))
            for w in ln:
                line_array.append(w)      
        #print len(line_array)
        
        for line in line_array:
            #print line
            strand=line.split("\t")[0]
            #print strand, options.strand
            if ( options.strand=='F' and strand=="F") or (options.strand=='R' and strand=="R"):
                #print "line %s" %line
                name=line.split("\t")[5]
                #print "name %s" %name
                if name!="NO_TU":
                    
                    start=(line.split("\t")[6])
                    #print start
                    end=(line.split("\t")[7])
                    Mvalue=float(line.split("\t")[8])
                    ttest=(line.split("\t")[9])
                    pvalue=(line.split("\t")[10].split("\n")[0])
                    if Mvalue>=options.thup or Mvalue<=options.thdown:
                        strand_array.append(strand)
                        nameinter_array.append(name)
                        #orfup_array.append(orfup)
                        #orfdown_array.append(orfdown)
                        #schema_array.append(schema)
                        start_array.append(int(start))
                        end_array.append(int(end))
                        #length_array.append(length)
                        Mvalue_array.append(float(Mvalue))
                        ttest_array.append(ttest)
                        pvalue_array.append(pvalue)
                        #if Mvalue>=0:
                        count_up=count_up+1
                    else:
                        count_down=count_down+1
        file_report.write("HD%s:\t%s\t%s\t\n" %(count_exp, count_up, count_down))
    
    
    
    
    matrice=[]
    for pippo in range(len(start_array)):
        array_interno=[]
        
        array_interno0=start_array[pippo]
        array_interno1=end_array[pippo]
        array_interno2= nameinter_array[pippo]
        array_interno3= strand_array[pippo]
        array_interno4=Mvalue_array[pippo]
        array_interno5=ttest_array[pippo]
        array_interno6=pvalue_array[pippo]
        
        array_interno.append(array_interno0)
        array_interno.append(array_interno1)
        array_interno.append(array_interno2)
        array_interno.append(array_interno3)
        array_interno.append(array_interno4)
        array_interno.append(array_interno5)
        array_interno.append(array_interno6)
        matrice.append(array_interno)
        
    #print matrice
    
    matrice.sort()
    i=0
    #orfup_array_sort=[]
    #orfdown_array_sort=[]
    nameinter_array_sort=[]
    strand_array_sort=[]
    #schema_array_sort=[]
    st_array=[]
    end_array=[]
    #length_array_sort=[]
    Mval_array=[]
    ttest_array_sort=[]
    pvalue_array_sort=[]
    for i in range(len(matrice)):
        
        st_array.append(int(matrice[i][0]))
        file_report.write("%s\t" %(matrice[i][0]))
        end_array.append(int(matrice[i][1]))
        file_report.write("%s\t" %(matrice[i][1]))
        nameinter_array_sort.append(matrice[i][2])
        strand_array_sort.append(matrice[i][3])
        Mval_array.append(int(matrice[i][4]))
        file_report.write("%s\n" %(matrice[i][4]))
        ttest_array_sort.append(float(matrice[i][5]))
        pvalue_array_sort.append(float(matrice[i][6]))
        
    #print len(orfup_array_sort)
    #print len(st_array)
    
    
    u=0
    
    
    count=0
    #for u in range(len(st_array)-1):
        
    print "\nVertex num"
    l1=len(st_array)
    l2=len(end_array)
    print l1
    g=igraph.Graph(l1)
    #print g.degree()
    #adj_matrix = make_matrix(l1, l2)
    i=0
    k=0
    for i in range(len(st_array)):
        st_ref=st_array[i]
        end_ref=end_array[i]
        controllo=0
        
        for k in range(len(st_array)):
            if k>i and controllo==0:
                st_mut=st_array[k]
                end_mut=end_array[k]
                if st_mut>=st_ref and st_mut<=end_ref or end_mut>=st_ref and end_mut<=end_ref:
            
                    g.add_edges([(i,k)])
                    #print i, k
                    #adj_matrix[i][k]=1
                else:
                    controllo==1
                    #adj_matrix[i][k]=0
    
    #print "\nmatrix"pytho
    #print adj_matrix
    
    #g=igraph.Graph.Adjacency(adj_matrix, mode="UNDIRECTED")
    
    #print "\nGraph:"
    #print g
    #igraph.plot(g)
    
    print ("\nIs connected: %s" %igraph.GraphBase.is_connected(g))
    
    i=0
    k=0
    j=0
    vert_max_btw=[]
    st_def=[]
    end_def=[]
    st_temp=[]
    end_temp=[]
    Mval_temp=[]
    if igraph.GraphBase.is_connected(g)==False or igraph.GraphBase.is_connected(g)==True:
        #print "entro"
        comp=igraph.Graph.components(g)
        
        #print igraph.Vertex.attributes(comp)
        print "num of components: %s" %(len(comp))
        
        for i in range(len(comp)):
            
            g_sub=igraph.clustering.VertexClustering.subgraph(comp, i)
            #print "\nsub graph %s: %s" %(i, g_sub)
            vertex_array=comp[i]
            #print "vertexes old graph: %s" %(comp[i])
            #print g_sub.vs
            #print "adj matrix: %s" %(igraph.Graph.get_adjacency(g_sub))
            betweenness_array=igraph.GraphBase.betweenness(g_sub)
            ln=len(betweenness_array)
            #print "betweenness for each vertex: %s" %(betweenness_array)
            max_btw=max(betweenness_array)
            for k in range(0,ln):
                if max_btw==betweenness_array[k]:
                    vert_max_btw.append(vertex_array[k])
                    #print vertex_array[k]
            #vertice_max_btw=betweenness_array.index(max_btw)
            #print "vertici max betweenness: %s" %(vertex_array[vertice_max_btw])
            #print len(vert_max_btw)
            #if len(vert_max_btw)==1:
                #print 'enro'
             #   st_def.append(st_array[vert_max_btw[0]])
                #print st_array[vert_max_btw[0]]
              #  end_def.append(end_array[vert_max_btw[0]])
                #print end_array[vert_max_btw[0]]
            
            
            #new method
            if len(vertex_array)>1:
                diff_array=[]
                provvisorio=[]
                #scommentare per riprendere la betweenness e commentare il ciclo dopo
                """
                for j in range(len(vert_max_btw)):
                    st_temp.append(st_array[vert_max_btw[j]])
                    #print end_array[vert_max_btw[j]]
                    end_temp.append(end_array[vert_max_btw[j]])
                    Mval_temp.append(Mval_array[vert_max_btw[j]])
                """
                for j in range(len(vertex_array)):
                    st_temp.append(st_array[vertex_array[j]])
                    #print end_array[vert_max_btw[j]]
                    end_temp.append(end_array[vertex_array[j]])
                    Mval_temp.append(abs(Mval_array[vertex_array[j]]))
                
                
                #dicotomization procedure
                """
                Mval_temp.sort()
                #print Mval_temp
                for elem in range(len(Mval_temp)-1):
                    diff=Mval_temp[elem]-Mval_temp[elem+1]
                    diff_array.append(abs(diff))
                #print diff_array
                    #file_out.write("\nelem %s, elem+1 %s, diff %s\n" %(dicotom[elem], dicotom[elem+1], diff))
                maxdiff=max(diff_array)
                #print maxdiff
                indice=diff_array.index(maxdiff)
                #file_out.write("max %s, indice %s\n" %(maxdiff, indice))
                uz=0
                for uz in range(indice+1, len(Mval_temp)):
                    provvisorio.append(Mval_temp[uz])
                    #file_out.write("provvisorio %s\n" %(dicotom[u]))
                """
                
                for elem in range(len(Mval_temp)):
                    if Mval_temp[elem]>=0.57 or Mval_temp[elem]<=-0.57:
                        provvisorio.append(Mval_temp[elem])
                
                if len(provvisorio)>0:    
                    #print provvisorio               
                    st_temp_dopo=[]
                    end_temp_dopo=[]
                    vz=0
                    zs=0
                    for vz in range(len(provvisorio)):
                        k_cont=0
                        for zs in range(len(Mval_array)):
                            if abs(Mval_array[zs])==provvisorio[vz] and k_cont==0:
                                st_temp_dopo.append(st_array[zs])
                                end_temp_dopo.append(end_array[zs])
                                k_cont=1
                    
                    
                    
                    
                    
                    
                    l1=len(st_temp_dopo)
                    l2=len(end_temp_dopo)
                    #print l1, l2
                    #print st_temp_dopo, end_temp_dopo
                    adj_matrix = make_matrix(l1, l2)
                    ii=0
                    kk=0
                    for ii in range(len(st_temp_dopo)):
                        st_ref=st_temp_dopo[ii]
                        end_ref=end_temp_dopo[ii]
                        for kk in range(len(st_temp_dopo)):
                            #if k>=i:
                                st_mut=st_temp_dopo[kk]
                                end_mut=end_temp_dopo[kk]
                                if st_mut>=st_ref and st_mut<=end_ref or end_mut>=st_ref and end_mut<=end_ref:
                                    adj_matrix[ii][kk]=1
                                else:
                                    adj_matrix[ii][kk]=0
                    
                    #print "\nmatrix"
                    #print adj_matrix
                    
                    gg=igraph.Graph.Adjacency(adj_matrix, mode="UNDIRECTED")
                    
                    print "\nGraph:"
                    print gg
                    #igraph.plot(g)
                    
                    print ("\nIs connected: %s" %igraph.GraphBase.is_connected(gg))
                    
                    if igraph.GraphBase.is_connected(gg)==True:
                        st_av=average_mod(st_temp_dopo, provvisorio)
                        end_av=average_mod(end_temp_dopo, provvisorio)
                        #print "st av %s, end av %s" %(st_av, end_av)
                        file_out.write("%s\t%s\t%s\t\n"%(strand_array_sort[u], st_av, end_av))
                        st_def.append(st_av)
                        end_def.append(end_av)
                    else:
                        compp=igraph.Graph.components(gg)
                    
                        #print igraph.Vertex.attributes(comp)
                        print "num of components: %s" %(len(compp))
                        iii=0
                        for iii in range(len(compp)):
                            g_subb=igraph.clustering.VertexClustering.subgraph(compp, iii)
                            vertex_array_iii=compp[iii]
                            jj=0
                            for jj in range(len(vertex_array_iii)):
                                st_temp_again.append(st_temp_dopo[vertex_array_iii[jj]])
                                #print end_array[vert_max_btw[j]]
                                end_temp_again.append(end_temp_dopo[vertex_array_iii[jj]])
                                Mval_temp_again.append(abs(provvisorio[vertex_array_iii[jj]]))
                    
                            st_av=average_mod(st_temp_again, Mval_temp_again)
                            end_av=average_mod(end_temp_again, Mval_temp_again)
                            file_out.write("%s\t%s\t%s\t\n"%(strand_array_sort[u], st_av, end_av))
                            #print "st av %s, end av %s" %(st_av, end_av)
                            st_def.append(st_av)
                            end_def.append(end_av)
                            st_temp_again=[]
                            end_temp_again=[]
                            Mval_temp_again=[]
                    
                    
                    
                    #st_av=average_mod(st_temp_dopo, provvisorio)
                    #end_av=average_mod(end_temp_dopo, provvisorio)
                    
                    
                       
                    #st_av=average_mod(st_temp, Mval_temp)
                    #end_av=average_mod(end_temp, Mval_temp)
                        #print "st av %s, end av %s" %(st_av, end_av)
                    #st_def.append(st_av)
                    #end_def.append(end_av)
            vert_max_btw=[]
            st_temp=[]
            end_temp=[]
            Mval_temp=[]
            st_temp_dopo=[]
            end_temp_dopo=[]
            Mval_temp_dopo=[]
            st_temp_again=[]
            end_temp_again=[]
            Mval_temp_again=[]
            
            
            
            
            
            
            
    #print "\n\nstart def: %s" %(st_def)
    #print "\n\nend def %s:" %(end_def)  
    
    st_array=[]
    end_array=[]
    Mval_array=[]
   
    
    print ("results stored into file '%s'. "
        % os.path.abspath(file_out_name) )
       


if __name__=='__main__':
    main()        
               
            
    
    
