# (a) Viterbi algorithm

seq="TTTAGCACCGGATGCGGTATCAATCCTGGTATCGTTAAAGCCTAGTGTTTCAAAAGTTCGAAAAACGGGCCGCGCAGCCCCGGGGCGTTCCCGGTAGGCTCGCGGCGGCTGGAGCACAGCGCCGCCGCGCCAGGGGACCCGCCGACTCCGTGGGCGTTCGGACCGCCTGTGTCCCTTCCGGAAGGCCTGGCAGGGTCGCCTCTTCCCGGAGCGGCCGCCCACACCGCGGCGCATGCAGGGCCCCGGTACGTTGACGTCTCTGACCCGCCGCTGCGGTGCGCCGGCGCGCCCGGCCGGGCGCCCCGGCGCGGGGCTCCCACGAGGCCCGCATCACGCGCGCCACCCAGGTTACCGGGCTCGGCCGGTGGACCATCGGGAGGGGCGGGGGCGGACGCCCGACCCCGCGGGCCACTCCAGTTTCTTAACTTGATAACGAAGCACTGAATACAGAAACAAGTTAAATCTCCCGGGTGTCGACGCGGTCCTGGTAACTTTTCAAGTAGGGTCGGTCAATGAAGGAGTGTTAGGCCTGGCCCGCCGGAGGGTTGAAGGACGGGTACGCGGTGTCAACGCCCGCCCCGGTCCGCCGAAGTCCCCCTTCACCGAGGGCCGAGCCGTCCTGCCCCCGGAGTTTCGGCGCGACCCCTACCGGGCTGGCCGCCGGGCGGGAGTCCCGGCTGTGGGCGGCTGGTCACCGGGTTTGGACCTCGGCAGGACAGCAGTTTGCATTCATATACAGGAAAAACACTCCCCACATGTGGTAAATACGTTGAACAAAGTATGTTACCATGAATAAGTAGCACAGCTAGTATTTTGTTGTTAGCAAGGAAAATAAATTCGTAAATTATTATTAATACGTATTTATTGAGTACTTTATTAAACTAAATATGTAGTTGATATCAGGCTGACACAACAACATGGTCTAGCAAAGTTCCTAACATTGGTATTTACCGCCACTTGCAGGCTCGAGGCGCCAGAGGGGCCCTGCGCTAGGTGCGGC"
seq=strsplit(seq,split = '')[[1]]

# emission probability
emission_h=c(0.13,0.37,0.37,0.13)
emission_l=c(0.32,0.18,0.18,0.32)


# transition probability
u=0.01

# states
states = c("h", "l")

# all the probability is log transformed

# transition matrixt
transition.matrix = (matrix(data = log(c(1-u, u, u, 1-u)), nrow = 2, ncol = 2, dimnames = list(c("h", "l"), c("h", "l"))))

# emission probabilities tell you what is the change of landing on each side given that the particular die is selected
emission.matrix = t(matrix(data = log(c(emission_h,emission_l)), nrow = 4, ncol = 2, dimnames = list(c('A','C','G','T'), c('h', 'l'))))

# the probability of begining with each state
a0 = log(c(0.5, 0.5))
names(a0)=c('h','l')

# store probabilities and track state
vki = data.frame()
state.track = NULL

# initialization:
prob=NULL
for (state in states){
    vk1 = a0[state] + emission.matrix[state,seq[1]]
    prob=c(prob,vk1)
}
vki = rbind(vki,prob)

state.track = rbind(state.track, as.character(states))

colnames(vki)  = states

# begin iteration
for (i in 2:length(seq)) {
    prob=NULL
    state.i_1=NULL
    for (state.l in states){# for l
        temp.prob=NULL
        state.pres=NULL
        for (state.k in states) {# for k
            prob.kl=vki[i-1,state.k] + transition.matrix[state.k,state.l]
            temp.prob=c(temp.prob,prob.kl)
        }
        
        max.k = max(temp.prob)
        # which state in i-1 ?
        state.pre = states[which(temp.prob == max(temp.prob))]
        state.i_1=c(as.character(state.i_1),as.character(state.pre))
        
        vk.prob = max.k + emission.matrix[state.l,seq[i]]
        prob=c(prob,vk.prob)
    }
    vki = rbind(vki,prob)
    state.track <- rbind(state.track, as.character(state.i_1))
    
}

colnames(state.track) = states

# find the path with highest probabilty

# highest probability at final position
maxfinal=max(vki[nrow(vki),])
# which state is it?
final.state=states[which(vki[nrow(vki),]==maxfinal)]
# trace back to find the best path
n = length(seq)
best.path=c(final.state)
for (i in n:2){
    # i-1 state
    final.state=state.track[i,final.state]
    best.path=c(final.state,best.path)
}

# print out the intervals
# find all the G-C rich positions
GC.index=which(best.path=='h')
start=GC.index[1]
end=NULL

GC.intervals=NULL
for (i in 2:length(GC.index)){
    if (GC.index[i]-GC.index[i-1]>1){
        GC.intervals=rbind(GC.intervals,c(start,GC.index[i-1]))
        start=GC.index[i]
    }
    
}
GC.intervals=rbind(GC.intervals,c(start,GC.index[length(GC.index)]))


# print out
linelength=50
block = floor(length(seq)/linelength)

for (l in 0:block){
    if (l < block){
        xsub = seq[seq((l*linelength+1),(l+1)*linelength,1)]
        ysub = best.path[seq((l*linelength+1),(l+1)*linelength,1)]
        
    }else{
        xsub = seq[seq((l*linelength+1),length(seq),1)]  
        ysub = best.path[seq((l*linelength+1),length(best.path),1)]  
        
    }
    xsub = paste(xsub,collapse = '')
    ysub = paste(ysub,collapse = '')
    cat(xsub,ysub,sep='\n')
    cat('\n')
}


# (b) Forward and backward algorithms
# define sum log probabilities function
sumLogProb=function(a,b){
    if (a > b){
        return(a + log1p(exp(b-a)))
    }else{
        return(b + log1p(exp(a-b))) 
    }
}

# forward probabilities
# initial:
fki = data.frame()

# initialization:
prob=NULL
for (state in states){
    fk1 = a0[state] + emission.matrix[state,seq[1]]
    prob=c(prob,fk1)
}
fki = rbind(fki,prob)
colnames(fki)  = states

# forward iteration
for (i in 2:length(seq)) {
    prob=NULL
    for (state.k in states){# for k
        # sum of fj(i-1)*ajk
        sum.temp=fki[i-1,]+transition.matrix[,state.k]
        sumlog = sumLogProb(sum.temp[1],sum.temp[2])
        prob.fki= emission.matrix[state.k,seq[i]] + sumlog
            
        prob=c(prob,prob.fki)
    }
    names(prob)=states
    fki = rbind(fki,prob)
}

# backward probabilities
# initial:
bki = data.frame()

# initialization:
prob=c(1,1)
bki = rbind(prob,bki)

colnames(bki)  = states

# backward iteration
for (i in (length(seq)-1):1) {
    prob=NULL
    for (state.k in states){# for k
        # sum of a*e*b
        sum.temp=bki[1,] + transition.matrix[state.k,] + emission.matrix[,seq[i+1]]
        sumlog = sumLogProb(sum.temp[1],sum.temp[2])
        prob.bki= sumlog
        
        prob=c(prob,prob.bki)
    }
    names(prob)=states
    bki = rbind(prob,bki)
}

# calculate h state probability
probh=NULL
fb=fki+bki
for (i in 1:length(seq)){
    temp=fb[i,]
    sumlog=sumLogProb(temp[1],temp[2])
    probh=c(probh,exp(fki[i,1]+bki[i,1]-sumlog))
}

# plot the result
par(mar=c(5.1, 4.1, 8.1, 4.1), xpd=TRUE)
plot(1:length(seq),probh,col="blue",pch=19,xlab="sequence position",ylab="P(x=h|y)")

# viterbi result
v.vec=rep(0.05,length(seq))
v.vec[best.path=='h']=0.95
points(1:length(seq),v.vec,col="green",pch=19)

legend("topright", inset=c(0,-0.5), legend=c("Posterior probability","Viterbi"), pch=c(19,19),col=c("blue","green"))
