gar.fun.neural<-function(out.var,mod.in,bar.plot=T,struct=NULL,x.lab=NULL,
                  y.lab=NULL, wts.only = F){
  
  require(ggplot2)
  require(plyr)
  
  # function works with neural networks from neuralnet, nnet, and RSNNS package
  # manual input vector of weights also okay
  
  #gets weights for neural network, output is list
  #if rescaled argument is true, weights are returned but rescaled based on abs value
  nnet.vals<-function(mod.in,nid,rel.rsc,struct.out=struct){
    
    require(scales)
    require(reshape)
    
    #neuralnet package
    if('nn' %in% class(mod.in)){
      struct.out<-unlist(lapply(mod.in$weights[[1]],ncol))
      struct.out<-struct.out[-length(struct.out)]
      struct.out<-c(
        length(mod.in$model.list$variables),
        struct.out,
        length(mod.in$model.list$response)
      )      	
      wts<-unlist(mod.in$weights[[1]])   
    }
    
    if(nid) wts<-rescale(abs(wts),c(1,rel.rsc))
    
    #convert wts to list with appropriate names 
    hid.struct<-struct.out[-c(length(struct.out))]
    row.nms<-NULL
    for(i in 1:length(hid.struct)){
      if(is.na(hid.struct[i+1])) break
      row.nms<-c(row.nms,rep(paste('hidden',i,seq(1:hid.struct[i+1])),each=1+hid.struct[i]))
    }
    row.nms<-c(
      row.nms,
      rep(paste('out',seq(1:struct.out[length(struct.out)])),each=1+struct.out[length(struct.out)-1])
    )
    out.ls<-data.frame(wts,row.nms)
    out.ls$row.nms<-factor(row.nms,levels=unique(row.nms),labels=unique(row.nms))
    out.ls<-split(out.ls$wts,f=out.ls$row.nms)
    
    assign('struct',struct.out,envir=environment(nnet.vals))
    
    out.ls
    
  }
  
  # get model weights
  best.wts<-nnet.vals(mod.in,nid=F,rel.rsc=5,struct.out=NULL)
  
  # weights only if T
  if(wts.only) return(best.wts)
  
  #get variable names from mod.in object
  #change to user input if supplied
  
  if('nn' %in% class(mod.in)){
    x.names<-mod.in$model.list$variables
    y.names<-mod.in$model.list$response
  }
  
  # get index value for response variable to measure
  if('numeric' %in% class(mod.in)){
    out.ind <-  as.numeric(gsub('^[A-Z]','',out.var))
  } else {
    out.ind<- grep(out.var, y.names)
  }
  
  #change variables names to user sub 
  if(!is.null(x.lab)){
    if(length(x.names) != length(x.lab)) stop('x.lab length not equal to number of input variables')
    else x.names<-x.lab
  }
  if(!is.null(y.lab)){
    y.names<-y.lab
  } else {
    y.names <- y.names[grep(out.var, y.names)]
  }
  
  # organize hidden layer weights for matrix mult
  inp.hid <- best.wts[grep('hidden', names(best.wts))]
  split_vals <- substr(names(inp.hid), 1, 8)
  inp.hid <- split(inp.hid, split_vals)
  inp.hid <- lapply(inp.hid, function(x) t(do.call('rbind', x))[-1, ])
  
  # final layer weights for output
  hid.out<-best.wts[[grep(paste('out',out.ind),names(best.wts))]][-1]
  
  # matrix multiplication of output layer with connecting hidden layer
  max_i <- length(inp.hid)
  sum_in <- as.matrix(inp.hid[[max_i]]) %*% matrix(hid.out)
  
  # recursive matrix multiplication for all remaining hidden layers
  # only for multiple hidden layers
  if(max_i != 1){ 
    
    for(i in (max_i - 1):1) sum_in <- as.matrix(inp.hid[[i]]) %*% sum_in
    
    # final contribution vector for all inputs
    inp.cont <- sum_in    
    
  } else {
    
    inp.cont <- sum_in
    
  }
  
  #get relative contribution
  #inp.cont/sum(inp.cont)
  rel.imp<-{
    signs<-sign(inp.cont)
    signs*rescale(abs(inp.cont),c(0,1))
  }
  
  if(!bar.plot){
    out <- data.frame(rel.imp)
    row.names(out) <- x.names
    return(out)
  }
  
  to_plo <- data.frame(rel.imp,x.names)[order(rel.imp),,drop = F]
  to_plo$x.names <- factor(x.names[order(rel.imp)], levels = x.names[order(rel.imp)])
  out_plo <- ggplot(to_plo, aes(x = x.names, y = rel.imp, fill = rel.imp,
                                colour = rel.imp)) + 
    geom_bar(stat = 'identity') + 
    scale_x_discrete(element_blank()) +
    scale_y_continuous(y.names)
  
  return(out_plo)
  
}

gar.fun.neural