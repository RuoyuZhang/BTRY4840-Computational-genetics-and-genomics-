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
transition.matrix = (matrix(data = log2(c(1-u, u, u, 1-u)), nrow = 2, ncol = 2, dimnames = list(c("h", "l"), c("h", "l"))))

# emission probabilities tell you what is the change of landing on each side given that the particular die is selected
emission.matrix = t(matrix(data = log2(c(emission_h,emission_l)), nrow = 4, ncol = 2, dimnames = list(c('A','C','G','T'), c('h', 'l'))))

# the probability of begining with each state
a0 = log2(c(0.5, 0.5))
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