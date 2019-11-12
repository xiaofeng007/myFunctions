



#DNA???򻥲?
revcom <- function(DNA){
  tmp <- chartr('acgtACGT','tgcaTGCA',DNA)
  tmp <- paste(rev(strsplit(tmp,"")[[1]]),collapse="")
  tmp
}


#??��???У??????????ݿ???ȡ???Ӻ???ǰ??N??
splitTail <- function(x,split="/",n=0){
	temp <- strsplit(as.character(x),split)
	sapply(temp,function(x) x[length(x)-n])
}

#ָʾ??????ֻ?ܽ???��ת??Ϊָʾ????
design <- function(dat=beefdat$breed){
  if(is.numeric(dat)) { levelV = dat } else { levelV = as.numeric(factor(dat)) }
  num <- length(levelV)
  X <- matrix(0,nrow = num,ncol=nlevels(factor(levelV)))
  for(i in 1:num){
    col <- levelV[i]
    X[i,col] <- 1
  }
  X
}

#??ȡASReml?????ļ?pvc?е?????
#dat?ǰ??ж?ȡ???ַ?????pattern??Ҫ??ȡ?ı?��????Ҫ??asreml??as?ļ??е?pin???ı?��??һ??
valueExtrPvc <- function(dat=pvc,pattern=" phenVar "){
  if(!is.character(dat)) stop("Input Data is not character, Please read Input File by readLines")
  temp <- strsplit(dat[grep(pattern,dat)]," ")[[1]]
  temp <- temp[temp!=""]
  temp[(length(temp)-1):length(temp)]
}

#??Ⱥֵ?±? #find outlier index
#dat??һ????��??nSd?Ǳ?׼?????????????ٸ???׼??????Ⱥֵ
#dat is a vector, nSd is the number of standard variance, which is the criterion to judge if it is a outlier
#???󷵻???Ⱥֵ???±?
#return the index of all outlier
idxOutlier <- function(dat,nSd=4){
  if(!is.vector(dat)) warning("dat is not vector, Please input a vector")
  dat <- as.numeric(dat)
  z=(dat-mean(dat,na.rm = T))/sd(dat,na.rm = T)
  idx <- which(z > nSd | z < -nSd )
  idx
}

#??ָ??????��???????ݵ?ĳһ?н???????
#??????һ??û???ظ??Ļ???ֱ????????????Ȼ??ֱ?ӵ??ü???
#???????ظ???????speSort
#dat?????ݣ?col ?????????????ţ?order ????Ҫ?ųɵ?˳????��
#dat refers to data; col refers to colnames or the order number; order refers to a specific vector whatever you want
speSort <- function(dat,col,order){
	res <- NULL
	for(i in 1:length(order)){
		temp <- dat[dat[,col]==order[i],]
		res <- rbind(res,temp)
	}
	res
}


#?ҳ?x??��????ͬ??��????һ????????
closedIdx_dis = function(x,leftPos,rightPos){
    if(length(leftPos) != length(rightPos)) stop("the length of leftPos and rightPos must be same!")

    # if x is within specific gene
    temp = which(leftPos<=x & x<=rightPos)
    if(length(temp) == 1) c(temp,0) else {    # if x is intergenic gene
        x_left = abs(x-leftPos)
        x_right = abs(x-rightPos)
        if(min(x_right) < min(x_left)) {idx = which.min(x_right);distance=min(x_right)} else{idx=which.min(x_left);distance=min(x_left)}
        c(idx,distance)
    }
}

findClosed = function(dat,ref){
    # ref must contains colummns below: "chr" for snp chromosome; "id" for gene around specific snp; "start" for start position of gene or QTL; "end" for end position of gene or QTL;
    # dat must contains columns below:  "chr" for snp chromosome; "bp" for position of snps
    if(!is.data.frame(dat)) stop("Input variable must be data.frame format! ")
    dat$id = NULL
    dat$dis = NULL
    for (i in 1:nrow(dat)){
        ref_temp = ref[ref$chr==dat[i,"chr"] ,]
        idx_dis = closedIdx_dis(as.numeric(dat[i,"bp"]),as.numeric(ref_temp$start),as.numeric(ref_temp$end))
        dat[i,"id"] = ref_temp[idx_dis[1],"id"]
        dat[i,"dis"] = idx_dis[2]
    }
    dat
}



