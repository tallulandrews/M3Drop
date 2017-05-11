#M3Drop_test<-function(x, size, mu){
#	prob <- -1;
#	out <- .C("test", as.integer(x), as.double(size), as.double(mu), as.double(prob));
#	return(out[[4]])
#}


# pnbinom<-function (q, size, prob, mu, lower.tail = TRUE, log.p = FALSE) 
#{
#    if (!missing(mu)) {
#        if (!missing(prob)) 
#            stop("'prob' and 'mu' both specified")
#        .Call(C_pnbinom_mu, q, size, mu, lower.tail, log.p)
#    }
#    else .Call(C_pnbinom, q, size, prob, lower.tail, log.p)
#}

