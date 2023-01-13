# Add 3prime 1kb of genes without overlapped regions of other genes
# Chul Lee (clee03@rockefeller.edu)

# program: extand_3prime
fNAME_fna    = "GCF_003957565.2_bTaeGut1.4.pri_genomic.fna"
fNAME_gtf    = "GCF_003957565.2_bTaeGut1.4.pri_genomic_ManuallyAddingGene_tRNArRNA.gtf"
oNAME_gtf    = "GCF_003957565.2_bTaeGut1.4.pri_genomic_ext1kb.gtf"
oNAME_filteredgtf = "GCF_003957565.2_bTaeGut1.4.pri_genomic_ext1kb_filtered.gtf"

# Library
import sys
import os

if len(sys.argv) == 5:
    try:
        fNAME_fna = sys.argv[1]
        fNAME_gtf = sys.argv[2]
        oNAME_gtf = sys.argv[3]
        oNAME_filteredgtf = sys.argv[4]
    except:
        pass
        
# function
def nScaff_iLen_mat(fNAME_fna):
    nScaff_iLen_dic = {}
    fpin= open(fNAME_fna,'r')
    for line in fpin:
        if line[0]==">":
            nScaff = line[1:].strip('\n').split(" ")[0]
            nScaff_iLen_dic.setdefault(nScaff,0)
        else:
            line = line.strip('\n')
            line = line.replace(' ','')
            nScaff_iLen_dic[nScaff]+=len(line)
    fpin.close()
    return nScaff_iLen_dic

def ext1kb_gtf(fNAME_fna,fNAME_gtf,oNAME_gtf):
    nScaff_iLen_dic = nScaff_iLen_mat(fNAME_fna)
    #print(nScaff_iLen_dic)
    
    iFlag = -1
    fpin = open(fNAME_gtf,'r')
    fpout = open(oNAME_gtf,'w')
    for line in fpin:
        if line[0] == "#":
            if not line.strip('\n')=="###":
                fpout.write(line)
        else:
            part        = line.strip('\n').split("\t")
            nScaf       = part[0]
            nSource     = part[1]
            nType       = part[2]
            nStart      = part[3]
            nEnd        = part[4]
            nStrand     = part[6]
            nInfo_list  = part[8].split(";")

            if nType == "gene":
                if not iFlag == 1:
                    iFlag = 1
                    
                else:
                    for nTranscriptID in nTranscriptID_list:
                        # Transcript extension
                        tmpline         = nTranscriptID_nInfo_dic[nTranscriptID]
                        tmp_part        = tmpline.strip('\n').split('\t')
                        tmp_nScaf       = tmp_part[0]
                        tmp_nSource     = tmp_part[1]
                        tmp_nType       = tmp_part[2]
                        tmp_nStart      = tmp_part[3]
                        tmp_nEnd        = tmp_part[4]
                        tmp_nStrand     = tmp_part[6]
                        tmp_nInfo_list  = tmp_part[8].split(";")
                            
                        if tmp_nStrand == "+":
                            if nScaff_iLen_dic[tmp_nScaf] >= (int(tmp_nEnd)+1000):  
                                tmp_nEnd = str(int(tmp_nEnd)+1000)
                            else:
                                tmp_nEnd = str(nScaff_iLen_dic[tmp_nScaf])
                        else: # nStrand == "-":
                            if 1 <= (int(tmp_nStart)-1000):  
                                tmp_nStart = str(int(tmp_nStart)-1000)
                            else:
                                tmp_nStart = "1"
                        tmpline = tmp_nScaf +"\t"+ tmp_nSource +"\t"+ tmp_nType +"\t"+ tmp_nStart +"\t"+ tmp_nEnd +"\t"+ '.' +"\t"+ tmp_nStrand +"\t"+ '.' +"\t"+ ';'.join(tmp_nInfo_list) +'\n'
                        fpout.write(tmpline)

                        # Exon extension
                        if len(nTranscriptID_ExonList_dic[nTranscriptID]) > 0:  
                            if tmp_nStrand == "+":
                                tmpline     = nTranscriptID_ExonList_dic[nTranscriptID][-1]
                                tmp_part        = tmpline.strip('\n').split('\t')
                                tmp_nScaf       = tmp_part[0]
                                tmp_nSource     = tmp_part[1]
                                tmp_nType       = tmp_part[2]
                                tmp_nStart      = tmp_part[3]
                                if nScaff_iLen_dic[tmp_nScaf] >= (int(tmp_part[4])+1000):  
                                    tmp_nEnd = str(int(tmp_part[4])+1000)
                                else:
                                    tmp_nEnd = str(nScaff_iLen_dic[tmp_nScaf])
                                tmp_nStrand     = tmp_part[6]
                                tmp_nInfo_list  = tmp_part[8].split(";")
                                tmpline = tmp_nScaf +"\t"+ tmp_nSource +"\t"+ tmp_nType +"\t"+ tmp_nStart +"\t"+ tmp_nEnd +"\t"+ '.' +"\t"+ tmp_nStrand +"\t"+ '.' +"\t"+ ';'.join(tmp_nInfo_list) +'\n'
                                nTranscriptID_ExonList_dic[nTranscriptID][-1] = tmpline
                                
                            else:
                                tmpline     = nTranscriptID_ExonList_dic[nTranscriptID][-1]
                                tmp_part        = tmpline.strip('\n').split('\t')
                                tmp_nScaf       = tmp_part[0]
                                tmp_nSource     = tmp_part[1]
                                tmp_nType       = tmp_part[2]
                                if 1 <= (int(tmp_part[3])-1000):  
                                    tmp_nStart = str(int(tmp_part[3])-1000)
                                else:
                                    tmp_nStart = "1"
                                tmp_nEnd        = tmp_part[4]
                                tmp_nStrand     = tmp_part[6]
                                tmp_nInfo_list  = tmp_part[8].split(";")
                                tmpline = tmp_nScaf +"\t"+ tmp_nSource +"\t"+ tmp_nType +"\t"+ tmp_nStart +"\t"+ tmp_nEnd +"\t"+ '.' +"\t"+ tmp_nStrand +"\t"+ '.' +"\t"+ ';'.join(tmp_nInfo_list) +'\n'
                                nTranscriptID_ExonList_dic[nTranscriptID][-1] = tmpline

                            for ExonLine in nTranscriptID_ExonList_dic[nTranscriptID]:
                                fpout.write(ExonLine)

                        # Others without extension
                        if len(nTranscriptID_OtherList_dic[nTranscriptID]) > 0:
                            for OtherLine in nTranscriptID_OtherList_dic[nTranscriptID]:
                                fpout.write(OtherLine)
                
                nTranscriptID_list = []
                nTranscriptID_nInfo_dic = {}
                nTranscriptID_ExonList_dic = {}
                nTranscriptID_OtherList_dic = {}
                
                if nStrand == "+":
                    if nScaff_iLen_dic[nScaf] >= (int(nEnd)+1000):  
                        nEnd = str(int(nEnd)+1000)
                    else:
                        nEnd = str(nScaff_iLen_dic[nScaf])
                else: # nStrand == "-":
                    if 1 <= (int(nStart)-1000):  
                        nStart = str(int(nStart)-1000)
                    else:
                        nStart = "1"
                tmpline = nScaf +"\t"+ nSource +"\t"+ nType +"\t"+ nStart +"\t"+ nEnd +"\t"+ '.' +"\t"+ nStrand +"\t"+ '.' +"\t"+ ';'.join(nInfo_list) +'\n'
                fpout.write(tmpline)    

            
            elif nType == "transcript":
                nTranscriptID = nInfo_list[1]
                nTranscriptID_list.append(nTranscriptID)
                nTranscriptID_nInfo_dic.setdefault(nTranscriptID,line)
                nTranscriptID_ExonList_dic.setdefault(nTranscriptID,[])
                nTranscriptID_OtherList_dic.setdefault(nTranscriptID,[])

            elif nType == "exon":
                nTranscriptID_ExonList_dic[nTranscriptID].append(line)

            else:
                nTranscriptID_OtherList_dic[nTranscriptID].append(line)
    fpin.close()

    
    for nTranscriptID in nTranscriptID_list:
        # Transcript extension
        tmpline  = nTranscriptID_nInfo_dic[nTranscriptID]
        part = tmpline.strip('\n').split('\t')
        nScaf       = part[0]
        nSource     = part[1]
        nType       = part[2]
        nStart      = part[3]
        nEnd        = part[4]
        nStrand     = part[6]
        nInfo_list  = part[8].split(";")
                            
        if nStrand == "+":
            if nScaff_iLen_dic[nScaf] >= (int(nEnd)+1000):  
                nEnd = str(int(nEnd)+1000)
            else:
                nEnd = str(nScaff_iLen_dic[nScaf])
        else: # nStrand == "-":
            if 1 <= (int(nStart)-1000):  
                nStart = str(int(nStart)-1000)
            else:
                nStart = "1"
        tmpline = nScaf +"\t"+ nSource +"\t"+ nType +"\t"+ nStart +"\t"+ nEnd +"\t"+ '.' +"\t"+ nStrand +"\t"+ '.' +"\t"+ ';'.join(nInfo_list) +'\n'
        fpout.write(tmpline)

        # Exon extension
        if len(nTranscriptID_ExonList_dic[nTranscriptID]) > 0:  
            if nStrand == "+":
                tmpline     = nTranscriptID_ExonList_dic[nTranscriptID][-1]
                part        = tmpline.strip('\n').split('\t')
                nScaf       = part[0]
                nSource     = part[1]
                nType       = part[2]
                nStart      = part[3]
                if nScaff_iLen_dic[nScaf] >= (int(part[4])+1000):  
                    nEnd = str(int(part[4])+1000)
                else:
                    nEnd = str(nScaff_iLen_dic[nScaf])
                nStrand     = part[6]
                nInfo_list  = part[8].split(";")
                tmpline = nScaf +"\t"+ nSource +"\t"+ nType +"\t"+ nStart +"\t"+ nEnd +"\t"+ '.' +"\t"+ nStrand +"\t"+ '.' +"\t"+ ';'.join(nInfo_list) +'\n'
                nTranscriptID_ExonList_dic[nTranscriptID][-1] = tmpline
               
            else:
                tmpline     = nTranscriptID_ExonList_dic[nTranscriptID][-1]
                part        = tmpline.strip('\n').split('\t')
                nScaf       = part[0]
                nSource     = part[1]
                nType       = part[2]
                if 1 <= (int(part[3])-1000):  
                    nStart = str(int(part[3])-1000)
                else:
                    nStart = "1"
                nEnd        = part[4]
                nStrand     = part[6]
                nInfo_list  = part[8].split(";")
                tmpline = nScaf +"\t"+ nSource +"\t"+ nType +"\t"+ nStart +"\t"+ nEnd +"\t"+ '.' +"\t"+ nStrand +"\t"+ '.' +"\t"+ ';'.join(nInfo_list) +'\n'
                nTranscriptID_ExonList_dic[nTranscriptID][-1] = tmpline

            for ExonLine in nTranscriptID_ExonList_dic[nTranscriptID]:
                fpout.write(ExonLine)
        # Others without extension
        if len(nTranscriptID_OtherList_dic[nTranscriptID]) > 0:
            for OtherLine in nTranscriptID_OtherList_dic[nTranscriptID]:
                fpout.write(OtherLine)
    fpout.write("###\n")
    fpout.close()
                
def Exclude_Overlapped_Ext(fNAME_gtf, oNAME_gtf, oNAME_filteredgtf):
    nScaff_nGeneID_nInfo_dic = {}
    iCNT_gene = 0
    fpin = open(fNAME_gtf,'r')
    for line in fpin:
        if not line[0]=="#":
            part        = line.strip('\n').split("\t")
            nScaf       = part[0]
            nSource     = part[1]
            nType       = part[2]
            nStart      = part[3]
            nEnd        = part[4]
            nStrand     = part[6]
            nInfo_list  = part[8].split(";")

            if nType == "gene":
                nGeneID = nInfo_list[0]
                nScaff_nGeneID_nInfo_dic.setdefault(nScaf, {})
                nScaff_nGeneID_nInfo_dic[nScaf].setdefault(nGeneID, [nStart, nEnd, nStrand, ""])
                nScaff_nGeneID_nInfo_dic[nScaf][nGeneID][3]+=line
            else:
                try:
                    nScaff_nGeneID_nInfo_dic[nScaf][nGeneID][3]+=line
                except:
                    print(nScaf)
                    print(nScaff_nGeneID_nInfo_dic[nScaf].keys())
                    print(nGeneID)
                    print(line)
                    sys.exit()
    fpin.close()

    iFlag = 0
    fpin = open(oNAME_gtf,'r')
    fpout = open(oNAME_filteredgtf,'w')
    
    for line in fpin:
        if line[0] == "#":
            if not line.strip('\n')=="###":
                fpout.write(line)
        else:
            part        = line.strip('\n').split("\t")
            nScaf       = part[0]
            nSource     = part[1]
            nType       = part[2]
            nStart      = part[3]
            nEnd        = part[4]
            nStrand     = part[6]
            nInfo       = part[8].strip()
            nGene = nInfo.split(";")[0]
            iStart = int(nStart)
            iEnd   = int(nEnd)

            if nStrand == "+":
                ext_iStart = int(nEnd)-1000
                ext_iEnd   = int(nEnd)
            else:
                ext_iStart = int(nStart)
                ext_iEnd   = int(nStart) + 1000

            if nType == "gene":
                iFlag = 0
                for que_nGeneID in nScaff_nGeneID_nInfo_dic[nScaf].keys():
                    pass_flag = 0
                    que_iStart  = int(nScaff_nGeneID_nInfo_dic[nScaf][que_nGeneID][0])
                    que_iEnd    = int(nScaff_nGeneID_nInfo_dic[nScaf][que_nGeneID][1])
                    que_nStrand = nScaff_nGeneID_nInfo_dic[nScaf][que_nGeneID][2]
                    que_nInfo   = nScaff_nGeneID_nInfo_dic[nScaf][que_nGeneID][3].split('\n')[0].split("\t")[8]
                    if nInfo in que_nInfo:
                        pass_flag = 1
                    if pass_flag == 0:                        
                        if que_nStrand == "+":
                            que_iEnd += 1000
                        else:
                            que_iStart -= 1000                            
                        if que_iStart <= ext_iStart <= que_iEnd:
                            iFlag = 0
                            break
                        else:
                            if que_iStart <= ext_iEnd <= que_iEnd:
                                iFlag = 0
                                break
                            else:
                                iFlag = 1
                if iFlag == 1:
                    fpout.write(line)
                else:
                    tmpline = nScaff_nGeneID_nInfo_dic[nScaf][nGene][3]
                    fpout.write(tmpline)
            else:
                if iFlag == 1:
                    fpout.write(line)
    fpin.close()
    fpout.write("###\n")
    fpout.close()
            

def main():
    ext1kb_gtf(fNAME_fna,fNAME_gtf,oNAME_gtf)
    Exclude_Overlapped_Ext(fNAME_gtf, oNAME_gtf, oNAME_filteredgtf)


# main
if __name__ == "__main__":
    main()
