# BTRY 4840 Problem set2 Ruoyu Zhang rz253
# Problem2
#seq1='GACTT'
#seq2='GGCAATC'

# sequences to compare
seq1="GGGTGGGAAAATAGACCAATAGGCAGAGAGAGTCAGTGCCTATCAGAAACCCAAGAGTCTTCTCTGTCTCCACATGCCCAGTTTCTATTGGTCTCCTTAAACCTGTCTTGTAACCTTGATA"
seq2="AAAGGGAAACATAGACAGGGGACACTCAAAGTTAGTGCCTGCTGGAAAGCAGACCTCTGTCTCCAAGCACCCAACTTCTACTTGTGAGCTGCCTTGTAACCTGGATA"

# score matrix
score_mat=matrix(c(91,-114,-31,-123,
                   -114,100,-215,-31,
                   -31,-125,100,-114,
                   -123,-31,-114,91),nrow = 4)

# assign row names and col names of the score matrix
colnames(score_mat)=c('A','C','G','T')
rownames(score_mat)=colnames(score_mat)

############
#(b)
# function for global alignment with linear gap penalties
global_l = function(score,seq1,seq2,d=100){
    # convert string into vector
    seq1=strsplit(seq1,split = '')[[1]]
    seq2=strsplit(seq2,split = '')[[1]]
    
    # length of seq1 and seq2
    m = length(seq1)
    n = length(seq2)
    
    # initialize F and T matrix
    Fmat = matrix(rep(0,(n+1)*(m+1)),ncol = m+1)
    Tmat = matrix(rep(0,(n+1)*(m+1)),ncol = m+1)
    
    # initialize first row and first column
    for (i in seq(2,n+1,1)){
        Fmat[i,1] = Fmat[i-1,1]-d
        Tmat[i,1] = 'u'
    }
    for (j in seq(2,m+1,1)){
        Fmat[1,j] = Fmat[1,j-1]-d
        Tmat[1,j] = 'l'
    }
    
    # fill in the two matrix
    for (i in seq(2,n+1,1)){
        for (j in seq(2,m+1,1)){
            sab = score[seq1[j-1],seq2[i-1]] # find value from scroe matrix
            vg = Fmat[i-1,j-1] + sab
            vu = Fmat[i-1,j] - d
            vl = Fmat[i,j-1] - d
            Fmat[i,j]=max(c(vg,vu,vl))
            Tmat[i,j]=c('g','u','l')[which.max(c(vg,vu,vl))[1]]
        }
    }
    
    # traceback
    x=NULL; y=NULL # initialize alignments
    i = n + 1
    j = m + 1
    align_score = 0
    while (i > 1 | j > 1){
        if (Tmat[i,j] == 'g'){
            x=c(seq1[j-1],x)
            y=c(seq2[i-1],y)
            
            align_score = align_score + (Fmat[i,j] - Fmat[i-1,j-1])
            i = i-1
            j = j-1
        } else if (Tmat[i,j]=='l'){
            x=c(seq1[j-1],x)
            y=c('-',y)
            align_score = align_score - d
            j = j-1
        } else if (Tmat[i,j]=='u'){
            x=c('-',x)
            y=c(seq2[i-1],y)
            align_score = align_score - d
            i = i-1
        }
    }
    
    # organize output
    linelength=40
    block = floor(max(m,n)/linelength)
    
    for (l in 0:block){
        if (l < block){
            xsub = x[seq((l*linelength+1),(l+1)*linelength,1)]
            ysub = y[seq((l*linelength+1),(l+1)*linelength,1)]
            
        }else{
            xsub = x[seq((l*linelength+1),length(x),1)]  
            ysub = y[seq((l*linelength+1),length(x),1)]  
            
        }
        xsub = paste(xsub,collapse = '')
        ysub = paste(ysub,collapse = '')
        cat(xsub,ysub,sep='\n')
        cat('\n')
    }
    return(align_score)
}

# alignment with linear gap penalties
global_l(score = score_mat,seq1 = seq1, seq2=seq2, d=100)


##########################
#(c)
# function for global alignment with affine gap penalties
global_a = function(score,seq1,seq2,d=430,e=30){
    # convert string into vector
    seq1=strsplit(seq1,split = '')[[1]]
    seq2=strsplit(seq2,split = '')[[1]]
    
    # length of seq1 and seq2
    m = length(seq1)
    n = length(seq2)
    
    # initialize M, I and Trace matrix
    Mmat = matrix(rep(0,(n+1)*(m+1)),ncol = m+1)
    Ix = matrix(rep(0,(n+1)*(m+1)),ncol = m+1)
    Iy = matrix(rep(0,(n+1)*(m+1)),ncol = m+1)
    
    MT = matrix(rep(0,(n+1)*(m+1)),ncol = m+1) # trace matrix for M
    XT = matrix(rep(0,(n+1)*(m+1)),ncol = m+1) # trace matrix for Ix
    YT = matrix(rep(0,(n+1)*(m+1)),ncol = m+1) # trace matrix for Iy
    
    # initialize first row and first column
    for (i in seq(2,n+1,1)){
        Ix[i,1] = Ix[i-1,1] - e
        Mmat[i,1] = -1/0
        Iy[i,1] = -1/0
        MT[i,1] = 'X'
        XT[i,1] = 'X'
        YT[i,1] = 'X'
    }
    for (j in seq(2,m+1,1)){
        Iy[1,j] = Iy[1,j-1] - e
        Mmat[1,j] = -1/0
        Ix[1,j] = -1/0
        MT[1,j] = 'Y'
        XT[1,j] = 'Y'
        YT[1,j] = 'Y'
    }
    
    # fill in the two matrix
    for (i in seq(2,n+1,1)){
        for (j in seq(2,m+1,1)){
            sab = score[seq1[j-1],seq2[i-1]] # find value from scroe matrix
            Mmat[i,j] = max(c(Mmat[i-1,j-1]+sab,Ix[i-1,j-1]+sab,Iy[i-1,j-1]+sab))
            Ix[i,j] = max(c(Mmat[i-1,j]-d,Ix[i-1,j]-e))
            Iy[i,j] = max(c(Mmat[i,j-1]-d,Iy[i,j-1]-e))
            
            #print(c('M','X','Y')[which.max(c(Mmat[i-1,j-1]+sab,Ix[i-1,j-1]+sab,Iy[i-1,j-1]+sab))[1]])
            MT[i,j] = c('M','X','Y')[which.max(c(Mmat[i-1,j-1]+sab,Ix[i-1,j-1]+sab,Iy[i-1,j-1]+sab))[1]]
            XT[i,j] = c('M','X')[which.max(c(Mmat[i-1,j]-d,Ix[i-1,j]-e))[1]]
            YT[i,j] = c('M','Y')[which.max(c(Mmat[i,j-1]-d,Iy[i,j-1]-e))[1]]
            
        }
    }
    
    # traceback
    x=NULL; y=NULL # initialize alignments
    i = n + 1
    j = m + 1
    align_score = 0
    # find which matrix to start with
    Mcurrent = c('M','X','Y')[which.max(c(Mmat[i,j],Ix[i,j],Iy[i,j]))[1]]
    
    while (i > 1 | j > 1){
        
        if (Mcurrent == 'M'){ # look up MT
            x=c(seq1[j-1],x)
            y=c(seq2[i-1],y)
            
            if ((i-1)!=1 & (j-1)!=1){
                align_score = align_score + (Mmat[i,j] - Mmat[i-1,j-1])
            }
            
            # find out which trace matrix to use for next step
            Mcurrent = MT[i,j]    
            
            i = i-1
            j = j-1
            
        } else if (Mcurrent == 'Y'){ # go left
            x=c(seq1[j-1],x)
            y=c('-',y)
            if ((i-1)!=1 & (j-1)!=1){
                align_score = align_score - (Iy[i,j] - Iy[i,j-1])
            }
            
            # find out which trace matrix to use for next step
            Mcurrent = YT[i,j]
            
            j = j-1
            
        } else if (Mcurrent == 'X'){# go up
            x=c('-',x)
            y=c(seq2[i-1],y)
            if ((i-1)!=1 & (j-1)!=1){
                align_score = align_score - (Ix[i,j] - Ix[i-1,j])
            }
            
            # find out which trace matrix to use for next step
            Mcurrent = XT[i,j]
            
            i = i-1
            
        }
    }
    
    # organize output
    linelength=40
    block = floor(max(m,n)/linelength)
    
    for (l in 0:block){
        if (l < block){
            xsub = x[seq((l*linelength+1),(l+1)*linelength,1)]
            ysub = y[seq((l*linelength+1),(l+1)*linelength,1)]
            
        }else{
            xsub = x[seq((l*linelength+1),length(x),1)]  
            ysub = y[seq((l*linelength+1),length(x),1)]  
            
        }
        xsub = paste(xsub,collapse = '')
        ysub = paste(ysub,collapse = '')
        cat(xsub,ysub,sep='\n')
        cat('\n')
    }
    return(align_score)
}

# alignment with affine gap penalties
global_a(score = score_mat,seq1 = seq1, seq2=seq2)

##############################
#(d)
# Making the gap opening very cheap, d=30, e=10
global_a(score=score_mat,seq1=seq1,seq2=seq2,d=30,e=10)

# Making the gap opening very expensive, d=1000, e=300
global_a(score=score_mat,seq1=seq1,seq2=seq2,d=1000,e=300)

# Making the gap-open/gap-extension small, eg, 1, d=430, e=430:
global_a(score=score_mat,seq1=seq1,seq2=seq2,d=430,e=430)

# Making the gap-open/gap-extension large, eg, 100, d=400, e=4
global_a(score=score_mat,seq1=seq1,seq2=seq2,d=400,e=4)
